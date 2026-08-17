"""
End-to-end checks that the wrong alignment format fails helpfully.

Mixing SAM and BAM up is the easy mistake to make with these commands, because
``samtools view`` is a pipe away from either. These tests drive the real
executable so that they cover the whole path a user hits: the guard, the
translation into a click usage error, the exit code, and the text that lands on
stderr.

The unit-level behaviour lives in tests/unit/test_formats.py; what is checked
here is that it is actually wired up.
"""

from pathlib import Path
import subprocess

import pytest

from tests.cli.conftest import CLIRunner, assert_exit_code


@pytest.fixture
def cli_runner():
    """Provide CLI runner instance."""
    return CLIRunner()


@pytest.fixture
def test_data_dir():
    """Get path to the existing test data directory."""
    return Path(__file__).parent.parent / 'data'


class TestBamSuppliedWhereSamExpected:
    """`filter` and `extract` read SAM text and must say so."""

    @pytest.mark.parametrize('command', ['filter', 'extract'])
    def test_bam_as_argument_is_rejected(self, cli_runner, test_data_dir, command):
        """
        A BAM passed as the input argument fails with samtools advice.

        Before the format guard this surfaced as a UnicodeDecodeError
        traceback from inside the parser.

        Parameters
        ----------
        cli_runner : CLIRunner
            Subprocess runner fixture.
        test_data_dir : pathlib.Path
            Directory holding the mock test dataset.
        command : str
            Sub-command under test.
        """
        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                command,
                '--ref-idx',
                str(test_data_dir / 'test.fna.fai'),
                str(test_data_dir / 'test.bam'),
            ]
        )

        assert exit_code != 0
        assert 'UnicodeDecodeError' not in stderr
        assert 'Traceback' not in stderr
        assert 'samtools view -h' in stderr
        assert f'teloclip {command}' in stderr

    def test_bam_on_stdin_is_rejected(self, test_data_dir):
        """
        A BAM piped in on stdin is rejected just as a file argument is.

        Driven with subprocess directly rather than through CLIRunner: the
        runner sends stdin as text, which would re-encode the BAM and stop it
        being a BAM by the time it arrived. The bytes have to reach the process
        intact for this to test anything.

        Parameters
        ----------
        test_data_dir : pathlib.Path
            Directory holding the mock test dataset.
        """
        bam_bytes = (test_data_dir / 'test.bam').read_bytes()

        result = subprocess.run(
            [
                'teloclip',
                'filter',
                '--ref-idx',
                str(test_data_dir / 'test.fna.fai'),
            ],
            input=bam_bytes,
            capture_output=True,
        )
        stderr = result.stderr.decode('utf-8', errors='replace')

        assert result.returncode != 0
        assert 'Traceback' not in stderr
        assert 'samtools view -h' in stderr

    def test_sam_still_works(self, cli_runner, test_data_dir):
        """
        The guard does not consume the stream it inspected.

        Peeking at the underlying buffer must leave the position untouched, or
        the first records would be silently dropped.

        Parameters
        ----------
        cli_runner : CLIRunner
            Subprocess runner fixture.
        test_data_dir : pathlib.Path
            Directory holding the mock test dataset.
        """
        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                '--ref-idx',
                str(test_data_dir / 'test.fna.fai'),
                '--min-anchor',
                '50',
                str(test_data_dir / 'test.sam'),
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)
        # The SAM header is written through, so the stream was read from byte 0.
        assert stdout.startswith('@')
        assert len(stdout.splitlines()) > 1


class TestSamSuppliedWhereBamExpected:
    """`extend` needs a sorted, indexed BAM."""

    def test_sam_is_rejected_with_a_conversion_recipe(
        self, cli_runner, test_data_dir, tmp_path
    ):
        """
        A SAM gets the convert-sort-index recipe, not 'create the .bai'.

        The previous message advised running `samtools index` on the SAM,
        which cannot work.

        Parameters
        ----------
        cli_runner : CLIRunner
            Subprocess runner fixture.
        test_data_dir : pathlib.Path
            Directory holding the mock test dataset.
        tmp_path : pathlib.Path
            Pytest temporary directory fixture.
        """
        sam_copy = tmp_path / 'alignments.sam'
        sam_copy.write_bytes((test_data_dir / 'test.sam').read_bytes())

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            ['extend', str(sam_copy), str(test_data_dir / 'test.fna')]
        )

        combined = stdout + stderr
        assert exit_code != 0
        assert 'samtools view -b' in combined
        assert 'samtools sort' in combined
        # The old, misleading advice must be gone.
        assert f'samtools index {sam_copy}' not in combined
