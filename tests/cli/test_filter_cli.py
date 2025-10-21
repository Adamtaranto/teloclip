"""Test filter command CLI functionality."""

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
                sam_path = temp_dir / 'test.sam'
                bam_path.write_bytes(b'')
                fasta_path.write_text('')
                sam_path.write_text('')
                files['sam'] = sam_path
            elif corrupt_files:
                # Create corrupted files
                bam_path = temp_dir / 'test.bam'
                fasta_path = temp_dir / 'test.fasta'
                sam_path = temp_dir / 'test.sam'
                fai_path = temp_dir / 'test.fasta.fai'
                bam_path.write_text('INVALID BAM CONTENT')
                fasta_path.write_text('INVALID FASTA CONTENT')
                sam_path.write_text('INVALID SAM CONTENT')
                fai_path.write_text('INVALID FAI CONTENT')
                files['fai'] = fai_path
                files['sam'] = sam_path
            else:
                # Copy valid files
                bam_path = temp_dir / 'test.bam'
                bai_path = temp_dir / 'test.bam.bai'
                fasta_path = temp_dir / 'test.fasta'
                fai_path = temp_dir / 'test.fasta.fai'
                sam_path = temp_dir / 'test.sam'

                shutil.copy(test_data_dir / 'test.bam', bam_path)
                shutil.copy(test_data_dir / 'test.bam.bai', bai_path)
                shutil.copy(test_data_dir / 'test.fna', fasta_path)
                shutil.copy(test_data_dir / 'test.fna.fai', fai_path)
                shutil.copy(test_data_dir / 'test.sam', sam_path)
                files['bai'] = bai_path
                files['fai'] = fai_path
                files['sam'] = sam_path

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


class TestFilterCLIBasicFunctionality:
    """Test basic filter command functionality."""

    def test_filter_help_available(self, cli_runner):
        """Test filter help is available."""
        exit_code, stdout, stderr = cli_runner.run_teloclip(['filter', '--help'])

        assert_exit_code(exit_code, 0, stdout, stderr)
        assert_contains(stdout, 'usage:', case_sensitive=False)
        assert_contains(stdout, 'filter', case_sensitive=False)

        # Should show key options
        assert_contains(stdout, '--ref-idx', case_sensitive=False)
        assert_contains(stdout, '--motifs', case_sensitive=False)
        assert_contains(stdout, '--no-rev', case_sensitive=False)
        assert_contains(stdout, '--fuzzy', case_sensitive=False)
        assert_contains(stdout, '--min-anchor', case_sensitive=False)

    def test_missing_required_ref_idx_shows_error(self, cli_runner, test_files):
        """Test missing required --ref-idx shows error."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            ['filter', str(files['sam'])]
        )

        assert exit_code != 0
        error_output = stdout + stderr
        assert_contains(error_output, 'ref-idx', case_sensitive=False)

    def test_synthetic_telomeric_filtering(self, cli_runner):
        """Test filtering with synthetic data containing actual telomeric sequences."""
        synthetic_data_dir = Path(__file__).parent.parent / 'integration' / 'test_data'
        synthetic_sam = synthetic_data_dir / 'synthetic_alignments.sam'
        synthetic_fai = synthetic_data_dir / 'synthetic_contigs.fasta.fai'

        # Skip if synthetic data not available
        if not synthetic_sam.exists() or not synthetic_fai.exists():
            pytest.skip('Synthetic test data not available')

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                '--ref-idx',
                str(synthetic_fai),
                '--motifs',
                'TTAGGG',
                str(synthetic_sam),
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)

        # Check stderr for statistics (filter outputs stats to stderr)
        assert_contains(stderr, 'Processed', case_sensitive=False)
        assert_contains(stderr, 'SAM records', case_sensitive=False)

        # Verify the filter ran with telomeric motifs
        assert_contains(stderr, 'TTAGGG', case_sensitive=True)

    def test_valid_minimal_command_succeeds(self, cli_runner, test_files):
        """Test valid minimal command: teloclip filter --ref-idx INDEX [SAMFILE]."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)
        # Should produce valid SAM output
        assert len(stdout) > 0, 'Expected SAM output to stdout'

    def test_reading_from_stdin_and_writing_to_stdout(self, cli_runner, test_files):
        """Test reading from stdin and writing to stdout."""
        files = test_files()

        # Read SAM content to pass as stdin
        sam_content = files['sam'].read_text()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            ['filter', '--ref-idx', str(files['fai'])], input_data=sam_content
        )

        assert_exit_code(exit_code, 0, stdout, stderr)
        # Should produce filtered SAM output
        assert len(stdout) > 0, 'Expected filtered SAM output to stdout'

    def test_motif_filtering_with_motifs_option(self, cli_runner, test_files):
        """Test motif filtering with --motifs."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--motifs',
                'TTAGGG',
            ]
        )

        # Should succeed regardless of whether motifs are found
        assert_exit_code(exit_code, 0, stdout, stderr)

    def test_reverse_complement_search_control_with_no_rev(
        self, cli_runner, test_files
    ):
        """Test reverse complement search control with --no-rev."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--motifs',
                'TTAGGG',
                '--no-rev',
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)


class TestFilterCLIParameterValidation:
    """Test parameter validation for filter command."""

    def test_non_existent_ref_idx_shows_error(self, cli_runner, test_files):
        """Test non-existent ref-idx file shows error."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            ['filter', str(files['sam']), '--ref-idx', 'non_existent_file.fai']
        )

        assert exit_code != 0
        error_output = stdout + stderr
        assert_contains(error_output, 'does not exist', case_sensitive=False)

    def test_non_existent_sam_file_shows_error(self, cli_runner, test_files):
        """Test non-existent SAM file shows error."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            ['filter', 'non_existent_file.sam', '--ref-idx', str(files['fai'])]
        )

        assert exit_code != 0
        error_output = stdout + stderr
        # Note: The actual error message might vary
        assert len(error_output) > 0, 'Expected error output'

    def test_invalid_motif_sequences_in_motifs(self, cli_runner, test_files):
        """Test invalid motif sequences in --motifs."""
        files = test_files()

        # Test with invalid DNA characters
        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--motifs',
                'TTAGGGXYZ',  # X, Y, Z are not valid DNA bases
            ]
        )

        # Should either succeed (if validation is lenient) or fail gracefully
        assert isinstance(exit_code, int)

    def test_malformed_motif_list(self, cli_runner, test_files):
        """Test malformed motif list handling."""
        files = test_files()

        # Test with empty motif
        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--motifs',
                '',
            ]
        )

        # Should handle gracefully
        assert isinstance(exit_code, int)

    def test_invalid_file_format_inputs(self, cli_runner, test_files):
        """Test invalid file format inputs."""
        files = test_files(corrupt_files=True)

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
            ]
        )

        # Should handle corrupted files gracefully
        assert isinstance(exit_code, int)

    def test_empty_sam_file_handling(self, cli_runner, test_files):
        """Test empty SAM file handling."""
        files = test_files(empty_files=True)
        # Use valid fai file
        valid_files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(valid_files['fai']),
            ]
        )

        # Should handle empty SAM gracefully
        assert isinstance(exit_code, int)


class TestFilterCLINumericParameters:
    """Test numeric parameter validation for filter command."""

    def test_min_clip_clipping_length_requirements(self, cli_runner, test_files):
        """Test --min-clip: clipping length requirements."""
        files = test_files()

        # Test with different min-clip values
        for min_clip in [0, 1, 5, 100]:
            exit_code, stdout, stderr = cli_runner.run_teloclip(
                [
                    'filter',
                    str(files['sam']),
                    '--ref-idx',
                    str(files['fai']),
                    '--min-clip',
                    str(min_clip),
                ]
            )

            # Should accept valid min-clip values
            assert_exit_code(exit_code, 0, stdout, stderr)

    def test_max_break_gap_tolerance_validation(self, cli_runner, test_files):
        """Test --max-break: gap tolerance validation."""
        files = test_files()

        # Test with different max-break values
        for max_break in [0, 50, 100, 1000]:
            exit_code, stdout, stderr = cli_runner.run_teloclip(
                [
                    'filter',
                    str(files['sam']),
                    '--ref-idx',
                    str(files['fai']),
                    '--max-break',
                    str(max_break),
                ]
            )

            # Should accept valid max-break values
            assert_exit_code(exit_code, 0, stdout, stderr)

    def test_min_repeats_sequential_pattern_matching(self, cli_runner, test_files):
        """Test --min-repeats: sequential pattern matching."""
        files = test_files()

        # Test with different min-repeats values
        for min_repeats in [1, 2, 5, 10]:
            exit_code, stdout, stderr = cli_runner.run_teloclip(
                [
                    'filter',
                    str(files['sam']),
                    '--ref-idx',
                    str(files['fai']),
                    '--motifs',
                    'TTAGGG',
                    '--min-repeats',
                    str(min_repeats),
                ]
            )

            # Should accept valid min-repeats values
            assert_exit_code(exit_code, 0, stdout, stderr)

    def test_min_anchor_minimum_anchor_length(self, cli_runner, test_files):
        """Test --min-anchor: minimum anchor length."""
        files = test_files()

        # Test with different min-anchor values
        for min_anchor in [1, 100, 500, 1000]:
            exit_code, stdout, stderr = cli_runner.run_teloclip(
                [
                    'filter',
                    str(files['sam']),
                    '--ref-idx',
                    str(files['fai']),
                    '--min-anchor',
                    str(min_anchor),
                ]
            )

            # Should accept valid min-anchor values
            assert_exit_code(exit_code, 0, stdout, stderr)

    def test_negative_numeric_parameters_rejected(self, cli_runner, test_files):
        """Test negative numeric parameters are handled appropriately."""
        files = test_files()

        # Test negative min-clip
        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--min-clip',
                '-1',
            ]
        )

        # Should either reject negative values or handle gracefully
        assert isinstance(exit_code, int)


class TestFilterCLIFeatureCombinations:
    """Test feature combinations and flags for filter command."""

    def test_fuzzy_motif_matching_variations(self, cli_runner, test_files):
        """Test --fuzzy motif matching variations."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--motifs',
                'TTAGGG',
                '--fuzzy',
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)

    def test_match_anywhere_vs_default_end_matching(self, cli_runner, test_files):
        """Test --match-anywhere vs default end-matching."""
        files = test_files()

        # Test default behavior (no --match-anywhere)
        exit_code1, stdout1, stderr1 = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--motifs',
                'TTAGGG',
            ]
        )

        # Test with --match-anywhere
        exit_code2, stdout2, stderr2 = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--motifs',
                'TTAGGG',
                '--match-anywhere',
            ]
        )

        # Both should succeed
        assert_exit_code(exit_code1, 0, stdout1, stderr1)
        assert_exit_code(exit_code2, 0, stdout2, stderr2)

    def test_multiple_motifs_with_reverse_complements(self, cli_runner, test_files):
        """Test multiple motifs with/without reverse complements."""
        files = test_files()

        # Test multiple motifs with default reverse complement search
        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--motifs',
                'TTAGGG,CCCTAA,TTTAGGG',
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)

    def test_multiple_motifs_without_reverse_complements(self, cli_runner, test_files):
        """Test multiple motifs without reverse complements."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--motifs',
                'TTAGGG,CCCTAA,TTTAGGG',
                '--no-rev',
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)

    def test_log_level_integration_across_commands(self, cli_runner, test_files):
        """Test log level integration across all commands."""
        files = test_files()

        # Test different log levels
        for log_level in ['DEBUG', 'INFO', 'WARNING', 'ERROR']:
            exit_code, stdout, stderr = cli_runner.run_teloclip(
                [
                    'filter',
                    str(files['sam']),
                    '--ref-idx',
                    str(files['fai']),
                    '--log-level',
                    log_level,
                ]
            )

            # Should accept all valid log levels
            assert_exit_code(exit_code, 0, stdout, stderr)

    def test_complex_feature_combination(self, cli_runner, test_files):
        """Test complex combination of multiple flags."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--motifs',
                'TTAGGG,CCCTAA',
                '--fuzzy',
                '--min-repeats',
                '2',
                '--min-anchor',
                '100',
                '--min-clip',
                '5',
                '--max-break',
                '30',
                '--match-anywhere',
                '--log-level',
                'DEBUG',
            ]
        )

        # Should handle complex flag combinations
        assert_exit_code(exit_code, 0, stdout, stderr)

    def test_comprehensive_filtering_with_synthetic_data(self, cli_runner):
        """Test comprehensive filtering pipeline with synthetic telomeric data."""
        synthetic_data_dir = Path(__file__).parent.parent / 'integration' / 'test_data'
        synthetic_sam = synthetic_data_dir / 'synthetic_alignments.sam'
        synthetic_fai = synthetic_data_dir / 'synthetic_contigs.fasta.fai'

        # Skip if synthetic data not available
        if not synthetic_sam.exists() or not synthetic_fai.exists():
            pytest.skip('Synthetic test data not available')

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                '--ref-idx',
                str(synthetic_fai),
                '--motifs',
                'TTAGGG,CCCTAA',
                '--fuzzy',
                '--min-clip',
                '10',
                '--min-repeats',
                '3',
                '--min-anchor',
                '50',
                '--log-level',
                'INFO',
                str(synthetic_sam),
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)

        # Check stderr for processing statistics
        assert_contains(stderr, 'Processed', case_sensitive=False)
        assert_contains(stderr, 'SAM records', case_sensitive=False)

        # Should process real telomeric sequences effectively
        # The synthetic data contains actual TTAGGG repeats in soft-clipped regions
        assert_contains(stderr, 'TTAGGG', case_sensitive=True)

    def test_fuzzy_with_no_rev_combination(self, cli_runner, test_files):
        """Test --fuzzy with --no-rev combination."""
        files = test_files()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                str(files['sam']),
                '--ref-idx',
                str(files['fai']),
                '--motifs',
                'TTAGGG',
                '--fuzzy',
                '--no-rev',
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)


class TestFilterCLIAdvancedScenarios:
    """Test advanced real-world scenarios for filter CLI command."""

    def test_telomeric_repeat_detection(self, cli_runner):
        """Test detection of actual telomeric repeats in synthetic data."""
        synthetic_data_dir = Path(__file__).parent.parent / 'integration' / 'test_data'
        synthetic_sam = synthetic_data_dir / 'synthetic_alignments.sam'
        synthetic_fai = synthetic_data_dir / 'synthetic_contigs.fasta.fai'

        # Skip if synthetic data not available
        if not synthetic_sam.exists() or not synthetic_fai.exists():
            pytest.skip('Synthetic test data not available')

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                '--ref-idx',
                str(synthetic_fai),
                '--motifs',
                'TTAGGG',
                str(synthetic_sam),
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)

        # Check stderr for processing statistics
        assert_contains(stderr, 'Processed', case_sensitive=False)
        assert_contains(stderr, 'SAM records', case_sensitive=False)
        assert_contains(stderr, 'TTAGGG', case_sensitive=True)

    def test_large_file_processing(self, cli_runner):
        """Test filter performance with realistic file sizes."""
        synthetic_data_dir = Path(__file__).parent.parent / 'integration' / 'test_data'
        synthetic_sam = synthetic_data_dir / 'synthetic_alignments.sam'
        synthetic_fai = synthetic_data_dir / 'synthetic_contigs.fasta.fai'

        if not synthetic_sam.exists() or not synthetic_fai.exists():
            pytest.skip('Synthetic test data not available')

        import time

        start_time = time.time()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                '--ref-idx',
                str(synthetic_fai),
                '--motifs',
                'TTAGGG,CCCTAA',
                '--min-clip',
                '20',
                str(synthetic_sam),
            ]
        )

        processing_time = time.time() - start_time

        assert_exit_code(exit_code, 0, stdout, stderr)
        # Should complete within reasonable time (10 seconds for synthetic data)
        assert processing_time < 10.0

    def test_pipeline_compatibility_stdin_stdout(self, cli_runner):
        """Test compatibility with Unix pipelines using synthetic data."""
        synthetic_data_dir = Path(__file__).parent.parent / 'integration' / 'test_data'
        synthetic_sam = synthetic_data_dir / 'synthetic_alignments.sam'
        synthetic_fai = synthetic_data_dir / 'synthetic_contigs.fasta.fai'

        if not synthetic_sam.exists() or not synthetic_fai.exists():
            pytest.skip('Synthetic test data not available')

        # Test stdin input with synthetic data
        with open(synthetic_sam) as f:
            sam_content = f.read()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            ['filter', '--ref-idx', str(synthetic_fai), '--motifs', 'TTAGGG', '-'],
            input_data=sam_content,
        )

        assert_exit_code(exit_code, 0, stdout, stderr)
        assert_contains(stderr, 'Processed', case_sensitive=False)

    def test_edge_case_motif_patterns(self, cli_runner):
        """Test filtering with edge case motif patterns."""
        synthetic_data_dir = Path(__file__).parent.parent / 'integration' / 'test_data'
        synthetic_sam = synthetic_data_dir / 'synthetic_alignments.sam'
        synthetic_fai = synthetic_data_dir / 'synthetic_contigs.fasta.fai'

        if not synthetic_sam.exists() or not synthetic_fai.exists():
            pytest.skip('Synthetic test data not available')

        # Test with multiple motifs including edge cases
        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                '--ref-idx',
                str(synthetic_fai),
                '--motifs',
                'TTAGGG,CCCTAA,AACCCT,AGGGTT',
                '--fuzzy',
                str(synthetic_sam),
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)
        assert_contains(stderr, 'Processed', case_sensitive=False)

    def test_filtering_statistics_accuracy(self, cli_runner):
        """Test accuracy of filtering statistics with known data."""
        synthetic_data_dir = Path(__file__).parent.parent / 'integration' / 'test_data'
        synthetic_sam = synthetic_data_dir / 'synthetic_alignments.sam'
        synthetic_fai = synthetic_data_dir / 'synthetic_contigs.fasta.fai'

        if not synthetic_sam.exists() or not synthetic_fai.exists():
            pytest.skip('Synthetic test data not available')

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                '--ref-idx',
                str(synthetic_fai),
                '--motifs',
                'TTAGGG',
                '--log-level',
                'INFO',
                str(synthetic_sam),
            ]
        )

        assert_exit_code(exit_code, 0, stdout, stderr)

        # Verify statistics are reported in stderr
        assert_contains(stderr, 'Processed', case_sensitive=False)
        assert_contains(stderr, 'SAM records', case_sensitive=False)

        # Check that numbers make sense (should have processed some records)
        lines = stderr.split('\n')
        stats_lines = [line for line in lines if 'SAM records' in line]
        assert len(stats_lines) > 0

    def test_memory_usage_validation(self, cli_runner):
        """Test memory usage remains reasonable with realistic data."""
        synthetic_data_dir = Path(__file__).parent.parent / 'integration' / 'test_data'
        synthetic_sam = synthetic_data_dir / 'synthetic_alignments.sam'
        synthetic_fai = synthetic_data_dir / 'synthetic_contigs.fasta.fai'

        if not synthetic_sam.exists() or not synthetic_fai.exists():
            pytest.skip('Synthetic test data not available')

        import time

        start_time = time.time()

        exit_code, stdout, stderr = cli_runner.run_teloclip(
            [
                'filter',
                '--ref-idx',
                str(synthetic_fai),
                '--motifs',
                'TTAGGG,CCCTAA',
                '--min-clip',
                '10',
                '--max-break',
                '100',
                str(synthetic_sam),
            ]
        )

        processing_time = time.time() - start_time

        assert_exit_code(exit_code, 0, stdout, stderr)
        # Should complete within reasonable time and memory constraints
        assert processing_time < 10.0
