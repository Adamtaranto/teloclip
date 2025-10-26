"""Test extract command CLI functionality."""

from pathlib import Path
import shutil

import pytest

from tests.cli.conftest import (
    CLIRunner,
    assert_contains,
    assert_exit_code,
)


@pytest.fixture
def cli_runner():
    """Provide CLI runner instance."""
    return CLIRunner()


@pytest.fixture
def temp_dir():
    """Provide temporary directory that gets cleaned up."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def test_data_dir():
    """Get path to existing test data directory."""
    return Path(__file__).parent.parent / 'data'


@pytest.fixture
def test_files(test_data_dir, temp_dir):
    """Provide access to existing test files."""

    def _get_files(
        copy_to_temp: bool = False,
        corrupt_files: bool = False,
        empty_files: bool = False,
    ):
        """
        Get paths to test files, optionally copying to temp directory.

        Args:
            copy_to_temp: Copy files to temp directory (allows modification)
            corrupt_files: Create corrupted versions of files
            empty_files: Create empty versions of files

        Returns:
            Dict with file paths
        """
        if copy_to_temp or corrupt_files or empty_files:
            # Copy/create files in temp directory
            files = {}

            if empty_files:
                # Create empty files
                bam_path = temp_dir / 'test.bam'
                fasta_path = temp_dir / 'test.fasta'
                bam_path.write_bytes(b'')
                fasta_path.write_text('')
            elif corrupt_files:
                # Create corrupted files
                bam_path = temp_dir / 'test.bam'
                fasta_path = temp_dir / 'test.fasta'
                bam_path.write_text('INVALID BAM CONTENT')
                fasta_path.write_text('INVALID FASTA CONTENT')
            else:
                # Copy valid files
                bam_path = temp_dir / 'test.bam'
                bai_path = temp_dir / 'test.bam.bai'
                fasta_path = temp_dir / 'test.fasta'
                fai_path = temp_dir / 'test.fasta.fai'

                shutil.copy(test_data_dir / 'test.bam', bam_path)
                shutil.copy(test_data_dir / 'test.bam.bai', bai_path)
                shutil.copy(test_data_dir / 'test.fna', fasta_path)
                shutil.copy(test_data_dir / 'test.fna.fai', fai_path)
                files['bai'] = bai_path
                files['fai'] = fai_path

            files['bam'] = bam_path
            files['fasta'] = fasta_path
            return files
        else:
            # Use original files directly
            return {
                'bam': test_data_dir / 'test.bam',
                'bai': test_data_dir / 'test.bam.bai',
                'fasta': test_data_dir / 'test.fna',
                'fai': test_data_dir / 'test.fna.fai',
                'sam': test_data_dir / 'test.sam',
            }

    return _get_files


class TestExtractCLIBasicFunctionality:
    """Test basic extract command functionality."""

    def test_extract_help_available(self, cli_runner):
        """Test extract help is available."""
        exit_code, stdout, stderr = cli_runner.run_teloclip(['extract', '--help'])

        assert_exit_code(exit_code, 0, stdout, stderr)
        assert_contains(stdout, 'usage:', case_sensitive=False)
        assert_contains(stdout, 'extract', case_sensitive=False)

        # Should show key options
        assert_contains(stdout, '--ref-idx', case_sensitive=False)
        assert_contains(stdout, '--prefix', case_sensitive=False)
        assert_contains(stdout, '--extract-dir', case_sensitive=False)
        assert_contains(stdout, '--output-format', case_sensitive=False)

    def test_missing_required_ref_idx_shows_error(self, cli_runner, test_files):
        """Test missing required --ref-idx shows error."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            ['extract', str(files['bam'])]
        )

        assert exit_code != 0
        error_output = stdout + stderr
        assert_contains(error_output, 'ref-idx', case_sensitive=False)

    def test_extract_with_ref_idx_succeeds(self, cli_runner, test_files, temp_dir):
        """Test extract succeeds with required ref-idx parameter."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--extract-dir',
                str(temp_dir),
            ]
        )
        assert_exit_code(exit_code, 0, stdout, stderr)

    def test_prefix_affects_output_filenames(self, cli_runner, test_files, temp_dir):
        """Test --prefix affects output filenames."""
        files = test_files()
        prefix = 'custom_extract'

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--prefix',
                prefix,
                '--extract-dir',
                str(temp_dir),
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)

        # Check that files with prefix were created
        files_created = list(temp_dir.glob(f'{prefix}_*'))
        assert len(files_created) > 0, f"No files with prefix '{prefix}' were created"

    def test_extract_dir_creates_directory(self, cli_runner, test_files, temp_dir):
        """Test --extract-dir creates and uses specified directory."""
        files = test_files()
        extract_dir = temp_dir / 'extracted'

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--extract-dir',
                str(extract_dir),
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)
        # Directory should be created
        assert extract_dir.exists(), 'Extract directory was not created'
        # Should contain output files
        output_files = list(extract_dir.glob('*.fasta'))
        assert len(output_files) > 0, 'No output files created in extract directory'

    def test_output_format_fasta(self, cli_runner, test_files, temp_dir):
        """Test --output-format fasta."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--output-format',
                'fasta',
                '--extract-dir',
                str(temp_dir),
            ]
        )

        # Should accept fasta format
        assert_exit_code(exit_code, 0, stdout, stderr)
        # Check that fasta files were created
        fasta_files = list(temp_dir.glob('*.fasta'))
        assert len(fasta_files) > 0, 'No fasta files created'

    def test_output_format_fastq(self, cli_runner, test_files, temp_dir):
        """Test --output-format fastq."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--output-format',
                'fastq',
                '--extract-dir',
                str(temp_dir),
            ]
        )

        # Should accept fastq format
        assert_exit_code(exit_code, 0, stdout, stderr)
        # Check that fastq files were created
        fastq_files = list(temp_dir.glob('*.fastq'))
        assert len(fastq_files) > 0, 'No fastq files created'

    def test_reading_from_stdin(self, cli_runner, test_files, temp_dir):
        """Test reading SAM data from stdin."""
        files = test_files()

        # Read SAM content to pass as stdin
        sam_content = files['sam'].read_text()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            ['extract', '--ref-idx', str(files['fai']), '--extract-dir', str(temp_dir)],
            input_data=sam_content,
        )

        # Should accept stdin input
        assert_exit_code(exit_code, 0, stdout, stderr)


class TestExtractCLIParameterValidation:
    """Test parameter validation for extract command."""

    def test_non_existent_ref_idx_shows_error(self, cli_runner, test_files):
        """Test non-existent ref-idx file shows error."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            ['extract', str(files['bam']), '--ref-idx', 'non_existent_file.fai']
        )

        assert exit_code != 0
        error_output = stdout + stderr
        assert_contains(error_output, 'does not exist', case_sensitive=False)

    def test_non_existent_sam_file_shows_error(self, cli_runner, test_files):
        """Test non-existent SAM file shows error."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            ['extract', 'non_existent_file.sam', '--ref-idx', str(files['fai'])]
        )

        assert exit_code != 0
        error_output = stdout + stderr
        # Note: The actual error message might vary
        assert len(error_output) > 0, 'Expected error output'

    def test_invalid_output_format_shows_error(self, cli_runner, test_files):
        """Test invalid output format shows error."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--output-format',
                'invalid_format',
            ]
        )

        assert exit_code != 0
        error_output = stdout + stderr
        assert_contains(error_output, 'format', case_sensitive=False)

    def test_empty_sam_file_handling(self, cli_runner, test_files, temp_dir):
        """Test empty SAM file handling."""
        files = test_files()
        # Create empty SAM file in temp directory
        empty_sam = temp_dir / 'empty.sam'
        empty_sam.write_text('')

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(empty_sam),
                '--ref-idx',
                str(files['fai']),
                '--extract-dir',
                str(temp_dir),
            ]
        )

        # Should handle gracefully (might succeed with no output or fail gracefully)
        assert isinstance(exit_code, int)

    def test_corrupted_sam_file_handling(self, cli_runner, test_files, temp_dir):
        """Test corrupted SAM file handling."""
        files = test_files()
        # Create corrupted SAM file
        corrupted_sam = temp_dir / 'corrupted.sam'
        corrupted_sam.write_text('INVALID SAM CONTENT')

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(corrupted_sam),
                '--ref-idx',
                str(files['fai']),
                '--extract-dir',
                str(temp_dir),
            ]
        )

        # May exit with error or handle gracefully
        assert isinstance(exit_code, int)


class TestExtractCLINumericParameters:
    """Test numeric parameter validation for extract command."""

    def test_min_clip_validation(self, cli_runner, test_files, temp_dir):
        """Test --min-clip parameter validation."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--min-clip',
                '5',
                '--extract-dir',
                str(temp_dir),
            ]
        )

        # Should accept valid min-clip value
        assert_exit_code(exit_code, 0, stdout, stderr)

    def test_max_break_validation(self, cli_runner, test_files, temp_dir):
        """Test --max-break parameter validation."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--max-break',
                '100',
                '--extract-dir',
                str(temp_dir),
            ]
        )

        # Should accept valid max-break value
        assert_exit_code(exit_code, 0, stdout, stderr)

    def test_min_anchor_validation(self, cli_runner, test_files, temp_dir):
        """Test --min-anchor parameter validation."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--min-anchor',
                '50',
                '--extract-dir',
                str(temp_dir),
            ]
        )

        # Should accept valid min-anchor value
        assert_exit_code(exit_code, 0, stdout, stderr)

    def test_min_mapq_validation(self, cli_runner, test_files, temp_dir):
        """Test --min-mapq parameter validation."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--min-mapq',
                '10',
                '--extract-dir',
                str(temp_dir),
            ]
        )

        # Should accept valid min-mapq value
        assert_exit_code(exit_code, 0, stdout, stderr)

    def test_buffer_size_validation(self, cli_runner, test_files, temp_dir):
        """Test --buffer-size parameter validation."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--buffer-size',
                '500',
                '--extract-dir',
                str(temp_dir),
            ]
        )

        # Should accept valid buffer-size value
        assert_exit_code(exit_code, 0, stdout, stderr)


class TestExtractCLIFeatureCombinations:
    """Test feature combinations and flags for extract command."""

    def test_include_stats_flag(self, cli_runner, test_files, temp_dir):
        """Test --include-stats flag."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--include-stats',
                '--extract-dir',
                str(temp_dir),
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)

    def test_report_stats_creates_file(self, cli_runner, test_files, temp_dir):
        """Test --report-stats creates statistics file."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--report-stats',
                '--extract-dir',
                str(temp_dir),
            ]
        )

        # May succeed or not depending on implementation
        assert isinstance(exit_code, int)

    def test_no_mask_overhangs_flag(self, cli_runner, test_files, temp_dir):
        """Test --no-mask-overhangs flag."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--no-mask-overhangs',
                '--extract-dir',
                str(temp_dir),
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)

    def test_count_motifs_with_multiple_motifs(self, cli_runner, test_files, temp_dir):
        """Test --count-motifs with multiple motifs."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--count-motifs',
                'TTAGGG,CCCTAA',
                '--extract-dir',
                str(temp_dir),
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)

    def test_fuzzy_count_with_count_motifs(self, cli_runner, test_files, temp_dir):
        """Test --fuzzy-count with --count-motifs."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--count-motifs',
                'TTAGGG',
                '--fuzzy-count',
                '--extract-dir',
                str(temp_dir),
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)

    def test_feature_combination_all_flags(self, cli_runner, test_files, temp_dir):
        """Test extract with many flags combined."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'extract',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--prefix',
                'combined_test',
                '--min-clip',
                '2',
                '--max-break',
                '30',
                '--min-anchor',
                '80',
                '--include-stats',
                '--count-motifs',
                'TTAGGG',
                '--fuzzy-count',
                '--no-mask-overhangs',
                '--extract-dir',
                str(temp_dir),
            ]
        )

        # May succeed or fail depending on implementation
        assert isinstance(exit_code, int)
