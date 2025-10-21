"""Fixtures and helpers for CLI tests."""

from pathlib import Path
import subprocess
import tempfile
from typing import List, Optional, Tuple

import pytest


class CLIRunner:
    """Helper class for running CLI commands and capturing output."""

    def __init__(self, timeout: int = 30):
        """Initialize CLI runner with timeout."""
        self.timeout = timeout

    def run_command(
        self,
        cmd: List[str],
        cwd: Optional[Path] = None,
        input_data: Optional[str] = None,
        check: bool = False,
    ) -> Tuple[int, str, str]:
        """
        Run a command and return (exit_code, stdout, stderr).

        Args:
            cmd: Command and arguments to run
            cwd: Working directory for command
            input_data: Data to pass to stdin
            check: If True, raise CalledProcessError on non-zero exit

        Returns:
            Tuple of (exit_code, stdout, stderr)
        """
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout,
                cwd=cwd,
                input=input_data,
                check=check,
            )
            return result.returncode, result.stdout, result.stderr
        except subprocess.TimeoutExpired:
            return -1, '', f'Command timed out after {self.timeout} seconds'
        except subprocess.CalledProcessError as e:
            return e.returncode, e.stdout, e.stderr
        except FileNotFoundError:
            return -1, '', 'Command not found'

    def run_teloclip(
        self,
        args: List[str],
        cwd: Optional[Path] = None,
        input_data: Optional[str] = None,
        check: bool = False,
    ) -> Tuple[int, str, str]:
        """
        Run teloclip command with given arguments.

        Args:
            args: Arguments to pass to teloclip
            cwd: Working directory
            input_data: Data to pass to stdin
            check: If True, raise CalledProcessError on non-zero exit

        Returns:
            Tuple of (exit_code, stdout, stderr)
        """
        cmd = ['teloclip'] + args
        return self.run_command(cmd, cwd=cwd, input_data=input_data, check=check)


@pytest.fixture
def cli_runner():
    """Provide CLI runner instance."""
    return CLIRunner()


@pytest.fixture
def temp_dir():
    """Provide temporary directory that gets cleaned up."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def test_data_dir():
    """Provide path to CLI test data directory."""
    return Path(__file__).parent / 'data'


@pytest.fixture
def minimal_sam_content():
    """Provide minimal valid SAM file content."""
    return """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
@SQ\tSN:chr2\tLN:2000
read1\t0\tchr1\t100\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t*
read2\t0\tchr2\t500\t60\t25M25S\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t*
"""


@pytest.fixture
def minimal_fasta_content():
    """Provide minimal valid FASTA file content."""
    return """>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chr2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
"""


@pytest.fixture
def create_test_files(temp_dir, minimal_sam_content, minimal_fasta_content):
    """Create minimal test files for CLI testing."""

    def _create_files(
        with_bam: bool = False,
        with_indices: bool = True,
        corrupt_sam: bool = False,
        empty_files: bool = False,
    ):
        """
        Create test files with various configurations.

        Args:
            with_bam: Also create BAM file
            with_indices: Create index files (.bai, .fai)
            corrupt_sam: Create corrupted SAM file
            empty_files: Create empty files

        Returns:
            Dict with file paths
        """
        files = {}

        # Create SAM file
        sam_path = temp_dir / 'test.sam'
        if empty_files:
            sam_path.write_text('')
        elif corrupt_sam:
            sam_path.write_text('INVALID SAM CONTENT')
        else:
            sam_path.write_text(minimal_sam_content)
        files['sam'] = sam_path

        # Create FASTA file
        fasta_path = temp_dir / 'test.fasta'
        if empty_files:
            fasta_path.write_text('')
        else:
            fasta_path.write_text(minimal_fasta_content)
        files['fasta'] = fasta_path

        # Create indices if requested
        if with_indices and not empty_files:
            # Create minimal FASTA index
            fai_path = temp_dir / 'test.fasta.fai'
            fai_content = 'chr1\t120\t6\t60\t61\nchr2\t120\t132\t60\t61\n'
            fai_path.write_text(fai_content)
            files['fai'] = fai_path

            # Create BAM and index if requested
            if with_bam:
                bam_path = temp_dir / 'test.bam'
                bai_path = temp_dir / 'test.bam.bai'
                # Create dummy BAM and BAI files (not valid, but present)
                bam_path.write_bytes(b'BAM\x01')  # Minimal BAM header
                bai_path.write_bytes(b'BAI\x01')  # Minimal BAI content
                files['bam'] = bam_path
                files['bai'] = bai_path

        return files

    return _create_files


def assert_exit_code(exit_code: int, expected: int, stdout: str, stderr: str):
    """Assert exit code matches expected, with helpful error message."""
    if exit_code != expected:
        msg = f'Expected exit code {expected}, got {exit_code}'
        if stdout:
            msg += f'\nSTDOUT:\n{stdout}'
        if stderr:
            msg += f'\nSTDERR:\n{stderr}'
        pytest.fail(msg)


def assert_contains(text: str, expected: str, case_sensitive: bool = True):
    """Assert that text contains expected string."""
    if not case_sensitive:
        text = text.lower()
        expected = expected.lower()

    if expected not in text:
        pytest.fail(f"Expected '{expected}' in text, but got:\n{text}")


def assert_file_exists(path: Path, should_exist: bool = True):
    """Assert file existence matches expectation."""
    if should_exist and not path.exists():
        pytest.fail(f"Expected file {path} to exist, but it doesn't")
    elif not should_exist and path.exists():
        pytest.fail(f'Expected file {path} to not exist, but it does')


def assert_file_not_empty(path: Path):
    """Assert file exists and has content."""
    assert_file_exists(path, should_exist=True)
    if path.stat().st_size == 0:
        pytest.fail(f"Expected file {path} to have content, but it's empty")
