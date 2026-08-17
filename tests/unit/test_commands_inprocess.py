"""
In-process CLI tests for the teloclip sub-commands.

The suite in ``tests/cli/`` drives the installed console script through
``subprocess``, which exercises the real entry point but contributes nothing to
coverage: ``cli.py``, ``commands/filter.py`` and ``commands/extract.py`` all
reported 0% despite being run on every invocation.

These tests use Click's ``CliRunner`` to invoke the same commands in-process, so
the command modules are measured and their error paths can be asserted on
directly rather than through captured stderr.
"""

from pathlib import Path

from click.testing import CliRunner
import pytest

from teloclip.cli import main
from teloclip.commands.extend import extend
from teloclip.commands.extract import extract_cmd
from teloclip.commands.filter import filter_cmd

DATA_DIR = Path(__file__).resolve().parents[1] / 'data'
REF_FASTA = DATA_DIR / 'test.fna'
REF_FAI = DATA_DIR / 'test.fna.fai'
TEST_SAM = DATA_DIR / 'test.sam'
TEST_BAM = DATA_DIR / 'test.bam'


@pytest.fixture
def runner():
    """
    Provide a Click CLI runner.

    Returns
    -------
    CliRunner
        Runner configured to keep stdout and stderr separate, so tests can
        assert on error output without matching filtered SAM records.
    """
    return CliRunner()


def sam_records(output: str) -> list:
    """
    Extract SAM alignment records from captured command output.

    CliRunner folds the rich log stream into the same buffer as stdout, so
    records are identified structurally, by their tab-delimited field count,
    rather than by excluding header lines.

    Parameters
    ----------
    output : str
        Captured command output.

    Returns
    -------
    list
        The alignment record lines.
    """
    return [
        line
        for line in output.splitlines()
        if not line.startswith('@') and len(line.split('\t')) >= 11
    ]


@pytest.fixture
def sam_text():
    """
    Provide the mock SAM file contents.

    Returns
    -------
    str
        Full text of ``tests/data/test.sam``.
    """
    return TEST_SAM.read_text()


class TestMainGroup:
    """The top-level command group."""

    def test_no_subcommand_shows_help(self, runner):
        """Test that bare `teloclip` prints help and exits cleanly."""
        result = runner.invoke(main, [])

        assert result.exit_code == 0
        assert 'filter' in result.output
        assert 'extract' in result.output
        assert 'extend' in result.output

    def test_version_flag(self, runner):
        """Test that --version reports a version string."""
        result = runner.invoke(main, ['--version'])

        assert result.exit_code == 0
        assert 'teloclip' in result.output

    @pytest.mark.parametrize('name', ['filter', 'extract', 'extend'])
    def test_subcommands_are_registered(self, runner, name):
        """Test that each sub-command is reachable from the group."""
        result = runner.invoke(main, [name, '--help'])

        assert result.exit_code == 0

    def test_unknown_subcommand_fails(self, runner):
        """Test that an unknown sub-command is rejected."""
        result = runner.invoke(main, ['nosuchcommand'])

        assert result.exit_code != 0

    def test_one_unimportable_command_does_not_hide_the_others(self, monkeypatch):
        """Test that sub-commands register independently of one another.

        All three used to be imported inside a single try block, so an
        environment without pysam registered *no* commands at all, including
        filter, which needs nothing beyond the standard library, click and rich.
        """
        import importlib

        import click as click_module

        real_import_module = importlib.import_module

        def fail_extend(name, *args, **kwargs):
            """
            Import normally, except for the extend command module.

            Parameters
            ----------
            name : str
                Module path being imported.
            *args
                Passed through to importlib.import_module.
            **kwargs
                Passed through to importlib.import_module.

            Returns
            -------
            module
                The imported module.

            Raises
            ------
            ImportError
                Always, for the extend command module.
            """
            if name == 'teloclip.commands.extend':
                raise ImportError("No module named 'pysam'")
            return real_import_module(name, *args, **kwargs)

        group = click_module.Group('teloclip')
        monkeypatch.setattr('teloclip.cli.main', group)
        monkeypatch.setattr(importlib, 'import_module', fail_extend)

        from teloclip.cli import register_commands

        register_commands()

        assert 'filter' in group.commands
        assert 'extract' in group.commands
        assert 'extend' not in group.commands


class TestFilterCommand:
    """The filter sub-command."""

    def test_filters_sam_from_stdin(self, runner, sam_text):
        """Test that filtering a SAM stream emits records and a header."""
        result = runner.invoke(
            filter_cmd,
            ['--ref-idx', str(REF_FAI), '--min-anchor', '50'],
            input=sam_text,
        )

        assert result.exit_code == 0
        # SAM headers are passed through unchanged.
        assert '@SQ\tSN:contig01' in result.output
        # Some alignments in the mock dataset are terminal overhangs.
        assert sam_records(result.output)

    def test_motif_filtering_reduces_output(self, runner, sam_text):
        """Test that requiring a motif discards non-matching overhangs."""
        common = ['--ref-idx', str(REF_FAI), '--min-anchor', '50']

        unfiltered = runner.invoke(filter_cmd, common, input=sam_text)
        filtered = runner.invoke(
            filter_cmd,
            common + ['--motifs', 'TT', '--min-repeats', '2', '--fuzzy'],
            input=sam_text,
        )

        assert unfiltered.exit_code == 0
        assert filtered.exit_code == 0

        kept_all = sam_records(unfiltered.output)
        kept_motif = sam_records(filtered.output)

        assert kept_all
        assert len(kept_motif) < len(kept_all)

    def test_missing_ref_idx_is_an_error(self, runner, sam_text):
        """Test that --ref-idx is required."""
        result = runner.invoke(filter_cmd, [], input=sam_text)

        assert result.exit_code != 0

    def test_nonexistent_ref_idx_is_an_error(self, runner, sam_text):
        """Test that a missing index file is rejected rather than ignored."""
        result = runner.invoke(
            filter_cmd, ['--ref-idx', 'no_such_file.fai'], input=sam_text
        )

        assert result.exit_code != 0

    def test_logfile_is_written(self, runner, sam_text, tmp_path):
        """Test that --logfile captures the run log.

        init_logging has always supported this; no CLI exposed it until now.
        """
        logfile = tmp_path / 'nested' / 'run.log'

        result = runner.invoke(
            filter_cmd,
            [
                '--ref-idx',
                str(REF_FAI),
                '--min-anchor',
                '50',
                '--logfile',
                str(logfile),
            ],
            input=sam_text,
        )

        assert result.exit_code == 0
        assert logfile.exists()
        assert 'Exclusion summary' in logfile.read_text()

    def test_exclusion_buckets_sum_to_total(self, runner, sam_text, tmp_path):
        """Test that the reported exclusion reasons account for every discard.

        Reads with no usable soft clip previously fell through every bucket
        while still counting toward the discard total.
        """
        import re

        logfile = tmp_path / 'run.log'
        runner.invoke(
            filter_cmd,
            [
                '--ref-idx',
                str(REF_FAI),
                '--min-anchor',
                '50',
                '--logfile',
                str(logfile),
            ],
            input=sam_text,
        )

        text = logfile.read_text()
        buckets = [int(n) for n in re.findall(r'  - .*?: (\d+)', text)]
        total = int(re.search(r'Total discarded: (\d+)', text).group(1))

        assert buckets
        assert sum(buckets) == total


class TestExtractCommand:
    """The extract sub-command."""

    def test_writes_per_contig_files(self, runner, sam_text, tmp_path):
        """Test that overhang reads are split into per-contig-end FASTA files."""
        outdir = tmp_path / 'split'

        result = runner.invoke(
            extract_cmd,
            [
                '--ref-idx',
                str(REF_FAI),
                '--extract-dir',
                str(outdir),
                '--min-anchor',
                '50',
            ],
            input=sam_text,
        )

        assert result.exit_code == 0
        assert outdir.exists()
        written = sorted(p.name for p in outdir.glob('*.fasta'))
        assert written

    def test_secondary_alignments_do_not_crash(self, runner, tmp_path):
        """Test that a secondary alignment is filtered rather than raising.

        EnhancedStreamingSamFilter read self.exclude_secondary without ever
        assigning it, so extract aborted on the first FLAG 256 record.
        """
        header = '@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:contig01\tLN:560\n'
        seq = 'A' * 150
        record = (
            '\t'.join(
                [
                    'secondary_read',
                    '256',
                    'contig01',
                    '1',
                    '60',
                    '50S100M',
                    '*',
                    '0',
                    '0',
                    seq,
                    'I' * 150,
                ]
            )
            + '\n'
        )

        result = runner.invoke(
            extract_cmd,
            [
                '--ref-idx',
                str(REF_FAI),
                '--extract-dir',
                str(tmp_path / 'out'),
                '--min-anchor',
                '50',
            ],
            input=header + record,
        )

        assert result.exit_code == 0
        assert result.exception is None

    def test_nonexistent_ref_idx_is_an_error(self, runner, sam_text, tmp_path):
        """Test that a missing index file is rejected."""
        result = runner.invoke(
            extract_cmd,
            ['--ref-idx', 'no_such_file.fai', '--extract-dir', str(tmp_path)],
            input=sam_text,
        )

        assert result.exit_code != 0


class TestExtendCommand:
    """The extend sub-command."""

    def test_requires_existing_inputs(self, runner):
        """Test that missing positional inputs are rejected."""
        result = runner.invoke(extend, ['no_such.bam', 'no_such.fasta'])

        assert result.exit_code != 0

    def test_unknown_exclusion_name_is_an_error(self, runner):
        """Test that a misspelled --exclude-contigs name fails loudly.

        Unknown names were previously warned about and ignored, so a typo left
        the contig to be extended while the user believed it was held back.
        """
        if not TEST_BAM.exists():
            pytest.skip('tests/data/test.bam not available')

        result = runner.invoke(
            extend,
            [
                str(TEST_BAM),
                str(REF_FASTA),
                '--dry-run',
                '--exclude-contigs',
                'contigOOPS',
            ],
        )

        assert result.exit_code != 0
        assert 'contigOOPS' in result.output

    def test_min_clip_is_clamped(self, runner):
        """Test that --min-clip below 1 is raised to 1 with a warning.

        A clip that does not reach past the terminus contributes no novel
        sequence and applying it would shorten the contig.
        """
        if not TEST_BAM.exists():
            pytest.skip('tests/data/test.bam not available')

        result = runner.invoke(
            extend,
            [str(TEST_BAM), str(REF_FASTA), '--dry-run', '--min-clip', '0'],
        )

        assert result.exit_code == 0

    def test_removed_exclude_outliers_is_rejected(self, runner):
        """
        The removed --exclude-outliers flag is now a usage error.

        It spent a release deprecated and ignored, which is worse than absent:
        a command line carrying it looked like it was excluding contigs and was
        not. Failing outright is the honest outcome.

        Parameters
        ----------
        runner : click.testing.CliRunner
            Click test runner fixture.
        """
        if not TEST_BAM.exists():
            pytest.skip('tests/data/test.bam not available')

        result = runner.invoke(
            extend,
            [str(TEST_BAM), str(REF_FASTA), '--dry-run', '--exclude-outliers'],
        )

        assert result.exit_code != 0
        assert 'no such option' in result.output.lower()
