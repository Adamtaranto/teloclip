"""
Unit tests for enhanced extract functionality.

Tests the new extract_io module with efficient FASTA/FASTQ writing,
statistics tracking, and motif integration.
"""

from io import StringIO
from pathlib import Path
import tempfile
from unittest.mock import patch

import pytest

from teloclip.extract_io import (
    EfficientSequenceWriter,
    ExtractionStats,
    MultiFileSequenceWriter,
)


class TestExtractionStats:
    """Test ExtractionStats functionality."""

    def test_stats_initialization(self):
        """Test stats object initialization."""
        stats = ExtractionStats()
        assert stats.total_alignments == 0
        assert stats.left_overhangs == 0
        assert stats.right_overhangs == 0
        assert len(stats.contigs_processed) == 0
        assert len(stats.motif_matches) == 0

    def test_record_alignment_left(self):
        """Test recording left overhang alignment."""
        stats = ExtractionStats()
        motif_counts = {'TTAGGG': 2, 'CCCTAA': 1}

        stats.record_alignment('contig1', True, motif_counts)

        assert stats.total_alignments == 1
        assert stats.left_overhangs == 1
        assert stats.right_overhangs == 0
        assert 'contig1' in stats.contigs_processed
        assert 'contig1' in stats.contigs_with_left
        assert 'contig1' not in stats.contigs_with_right
        assert stats.motif_matches['TTAGGG'] == 2
        assert stats.motif_matches['CCCTAA'] == 1

    def test_record_alignment_right(self):
        """Test recording right overhang alignment."""
        stats = ExtractionStats()

        stats.record_alignment('contig2', False)

        assert stats.total_alignments == 1
        assert stats.left_overhangs == 0
        assert stats.right_overhangs == 1
        assert 'contig2' in stats.contigs_processed
        assert 'contig2' not in stats.contigs_with_left
        assert 'contig2' in stats.contigs_with_right

    def test_record_filter(self):
        """Test recording filtered alignments."""
        stats = ExtractionStats()

        stats.record_filter('quality')
        stats.record_filter('anchor')
        stats.record_filter('quality')

        assert stats.filter_counts['quality'] == 2
        assert stats.filter_counts['anchor'] == 1
        assert stats.total_filtered == 3

    def test_generate_report(self):
        """Test report generation."""
        stats = ExtractionStats()

        # Add some data
        stats.record_alignment('contig1', True, {'TTAGGG': 3})
        stats.record_alignment('contig1', False, {'CCCTAA': 1})
        stats.record_alignment('contig2', True)
        stats.record_filter('quality')
        stats.record_filter('anchor')

        report = stats.generate_report()

        assert 'Extraction Statistics Report' in report
        assert 'Total alignments processed: 3' in report
        assert 'Left overhangs: 2' in report
        assert 'Right overhangs: 1' in report


class TestEfficientSequenceWriter:
    """Test EfficientSequenceWriter functionality."""

    def test_writer_initialization(self):
        """Test writer initialization."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
            writer = EfficientSequenceWriter(tmp.name, 'fasta', 100)
            assert writer.output_path == Path(tmp.name)
            assert writer.output_format == 'fasta'
            assert writer.buffer_size == 100
            assert len(writer.buffer) == 0
            assert writer.sequences_written == 0

    def test_writer_invalid_format(self):
        """Test writer with invalid format."""
        with pytest.raises(ValueError, match='Unsupported output format'):
            EfficientSequenceWriter(output_format='invalid')

    def test_write_fasta_to_file(self):
        """Test writing FASTA sequences to file."""
        with tempfile.NamedTemporaryFile(
            mode='w', delete=False, suffix='.fasta'
        ) as tmp:
            tmp_path = tmp.name

        with EfficientSequenceWriter(tmp_path, 'fasta', buffer_size=2) as writer:
            # Write some sequences
            writer.write_sequence('seq1', 'ATCG', 'test sequence 1')
            writer.write_sequence('seq2', 'GCTA', 'test sequence 2')
            # Should flush automatically when buffer fills
            writer.write_sequence('seq3', 'TTTT', 'test sequence 3')

        # Read back and verify
        with open(tmp_path, 'r') as f:
            content = f.read()

        assert '>seq1 test sequence 1' in content
        assert '>seq2 test sequence 2' in content
        assert '>seq3 test sequence 3' in content
        assert 'ATCG' in content
        assert 'GCTA' in content
        assert 'TTTT' in content

    def test_write_with_stats(self):
        """Test writing sequences with statistics in headers."""
        with tempfile.NamedTemporaryFile(
            mode='w', delete=False, suffix='.fasta'
        ) as tmp:
            tmp_path = tmp.name

        stats = {
            'mapq': 30,
            'clip_length': 25,
            'overhang_length': 25,
            'motif_counts': {'TTAGGG': 2, 'CCCTAA': 1},
        }

        with EfficientSequenceWriter(tmp_path, 'fasta', buffer_size=1) as writer:
            writer.write_sequence(
                'test_read', 'ATCGATCG', 'overhang sequence', stats=stats
            )

        with open(tmp_path, 'r') as f:
            content = f.read()

        assert 'mapq=30' in content
        assert 'clip_len=25' in content
        assert 'overhang_len=25' in content
        assert 'motifs=TTAGGG:2,CCCTAA:1' in content

    @patch('sys.stdout', new_callable=StringIO)
    def test_write_to_stdout(self, mock_stdout):
        """Test writing to stdout."""
        with EfficientSequenceWriter(
            output_path=None, output_format='fasta', buffer_size=1
        ) as writer:
            writer.write_sequence('seq1', 'ATCG', 'test sequence')

        output = mock_stdout.getvalue()
        assert '>seq1 test sequence' in output
        assert 'ATCG' in output


class TestMultiFileSequenceWriter:
    """Test MultiFileSequenceWriter functionality."""

    def test_multi_writer_initialization(self):
        """Test multi-file writer initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            writer = MultiFileSequenceWriter(tmpdir, 'test', 'fasta', 100)
            assert writer.base_dir == Path(tmpdir)
            assert writer.prefix == 'test'
            assert writer.output_format == 'fasta'
            assert writer.buffer_size == 100
            assert len(writer.file_handles) == 0

    def test_get_writer_creates_files(self):
        """Test that _get_file_handle creates appropriate file handles."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with MultiFileSequenceWriter(tmpdir, 'sample', 'fasta', 100) as writer:
                # Get file handles for different contigs/ends
                handle1 = writer._get_file_handle('contig1', 'L')
                handle2 = writer._get_file_handle('contig1', 'R')
                handle3 = writer._get_file_handle('contig2', 'L')

                # Should create separate file handles
                assert len(writer.file_handles) == 3
                assert ('contig1', 'L') in writer.file_handles
                assert ('contig1', 'R') in writer.file_handles
                assert ('contig2', 'L') in writer.file_handles

                # Verify file handles are properly stored
                assert writer.file_handles[('contig1', 'L')] is handle1
                assert writer.file_handles[('contig1', 'R')] is handle2
                assert writer.file_handles[('contig2', 'L')] is handle3

                # Should reuse existing file handles
                handle1_again = writer._get_file_handle('contig1', 'L')
                handle2_again = writer._get_file_handle('contig1', 'R')
                assert handle1 is handle1_again
                assert handle2 is handle2_again
                assert len(writer.file_handles) == 3

    def test_write_sequences_to_separate_files(self):
        """Test writing sequences to separate files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with MultiFileSequenceWriter(tmpdir, 'test', 'fasta', 1) as writer:
                # Write sequences for different contigs/ends
                writer.write_sequence(
                    'contig1', 'L', 'read1', 'ATCGATCG', 'left overhang'
                )
                writer.write_sequence(
                    'contig1', 'R', 'read2', 'GCTAGCTA', 'right overhang'
                )
                writer.write_sequence(
                    'contig2', 'L', 'read3', 'TTTTAAAA', 'another left'
                )

            # Check that files were created
            expected_files = [
                Path(tmpdir) / 'test_contig1_L.fasta',
                Path(tmpdir) / 'test_contig1_R.fasta',
                Path(tmpdir) / 'test_contig2_L.fasta',
            ]

            for file_path in expected_files:
                assert file_path.exists()

            # Check file contents
            with open(expected_files[0], 'r') as f:
                content = f.read()
                assert '>read1 left overhang' in content
                assert 'ATCGATCG' in content

            with open(expected_files[1], 'r') as f:
                content = f.read()
                assert '>read2 right overhang' in content
                assert 'GCTAGCTA' in content

    def test_no_prefix(self):
        """Test multi-writer without prefix."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with MultiFileSequenceWriter(tmpdir, None, 'fasta', 1) as writer:
                writer.write_sequence('contig1', 'L', 'read1', 'ATCG', 'test')

            # Should create file without prefix
            expected_file = Path(tmpdir) / 'contig1_L.fasta'
            assert expected_file.exists()


class TestIntegration:
    """Integration tests for extract functionality."""

    def test_stats_and_writer_integration(self):
        """Test integration between stats and writer."""
        stats = ExtractionStats()

        with tempfile.TemporaryDirectory() as tmpdir:
            with MultiFileSequenceWriter(tmpdir, 'sample', 'fasta', 1) as writer:
                # Simulate processing some alignments
                alignments = [
                    {'contig_name': 'chr1', 'end': 'L', 'motif_counts': {'TTAGGG': 1}},
                    {'contig_name': 'chr1', 'end': 'R', 'motif_counts': {'CCCTAA': 2}},
                    {'contig_name': 'chr2', 'end': 'L', 'motif_counts': None},
                ]

                for i, alignment in enumerate(alignments):
                    # Write sequence
                    writer.write_sequence(
                        alignment['contig_name'],
                        alignment['end'],
                        f'read_{i}',
                        'ATCGATCG',
                        f'test sequence {i}',
                    )

                    # Record stats
                    stats.record_alignment(
                        alignment['contig_name'],
                        alignment['end'] == 'L',
                        alignment['motif_counts'],
                    )

                # Files are automatically written in the new implementation

            # Verify stats
            assert stats.total_alignments == 3
            assert stats.left_overhangs == 2
            assert stats.right_overhangs == 1
            assert len(stats.contigs_processed) == 2
            assert stats.motif_matches['TTAGGG'] == 1
            assert stats.motif_matches['CCCTAA'] == 2

            # Verify files created
            expected_files = [
                Path(tmpdir) / 'sample_chr1_L.fasta',
                Path(tmpdir) / 'sample_chr1_R.fasta',
                Path(tmpdir) / 'sample_chr2_L.fasta',
            ]

            for file_path in expected_files:
                assert file_path.exists(), (
                    f'Expected file {file_path} was not created'
                )  # Generate and check report
        # Generate report with reference contigs to see contig analysis
        reference_contigs = {'chr1', 'chr2', 'chr3'}
        report = stats.generate_report(reference_contigs)
        assert 'Total alignments processed: 3' in report
        assert 'Contigs with overhangs: 2' in report
