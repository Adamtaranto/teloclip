"""Unit tests for teloclip.io.streaming module.

Tests file validation and streaming utilities.
"""

import gzip
from pathlib import Path
import tempfile

from teloclip.io.streaming import validate_indexed_files


def write_stub_bam(path: Path) -> None:
    """
    Write a file that passes the BAM format check.

    Validation now inspects the leading bytes before looking for the index, so
    a zero-byte file created with touch() is rejected as not being a BAM. Only
    the BGZF (gzip) magic number is checked, so an empty gzip member is enough
    to stand in for a real BAM here without pulling in pysam.

    Parameters
    ----------
    path : Path
        Where to write the stub.
    """
    path.write_bytes(gzip.compress(b''))


class TestValidateIndexedFiles:
    """Test file validation functions."""

    def test_validate_indexed_files_all_exist(self):
        """Test validation when all files and indexes exist."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create test files
            fasta_path = temp_path / 'test.fasta'
            bam_path = temp_path / 'test.bam'
            fai_path = temp_path / 'test.fasta.fai'
            bai_path = temp_path / 'test.bam.bai'

            for file_path in [fasta_path, fai_path, bai_path]:
                file_path.touch()
            write_stub_bam(bam_path)

            is_valid, error_msg = validate_indexed_files(fasta_path, bam_path)

            assert is_valid is True
            assert error_msg == ''

    def test_validate_indexed_files_fasta_missing(self):
        """Test validation when FASTA file is missing."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            fasta_path = temp_path / 'missing.fasta'
            bam_path = temp_path / 'test.bam'
            bam_path.touch()

            is_valid, error_msg = validate_indexed_files(fasta_path, bam_path)

            assert is_valid is False
            assert f'FASTA file not found: {fasta_path}' in error_msg

    def test_validate_indexed_files_bam_missing(self):
        """Test validation when BAM file is missing."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            fasta_path = temp_path / 'test.fasta'
            bam_path = temp_path / 'missing.bam'
            fasta_path.touch()

            is_valid, error_msg = validate_indexed_files(fasta_path, bam_path)

            assert is_valid is False
            assert f'BAM file not found: {bam_path}' in error_msg

    def test_validate_indexed_files_fai_missing(self):
        """Test validation when FASTA index is missing."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            fasta_path = temp_path / 'test.fasta'
            bam_path = temp_path / 'test.bam'
            bai_path = temp_path / 'test.bam.bai'

            # Create files but not .fai index
            for file_path in [fasta_path, bai_path]:
                file_path.touch()
            write_stub_bam(bam_path)

            is_valid, error_msg = validate_indexed_files(fasta_path, bam_path)
            expected_fai = temp_path / 'test.fasta.fai'

            assert is_valid is False
            assert f'FASTA index not found: {expected_fai}' in error_msg
            assert 'samtools faidx' in error_msg

    def test_validate_indexed_files_bai_missing(self):
        """Test validation when BAM index is missing."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            fasta_path = temp_path / 'test.fasta'
            bam_path = temp_path / 'test.bam'
            fai_path = temp_path / 'test.fasta.fai'

            # Create files but not .bai index
            for file_path in [fasta_path, fai_path]:
                file_path.touch()
            write_stub_bam(bam_path)

            is_valid, error_msg = validate_indexed_files(fasta_path, bam_path)
            expected_bai = temp_path / 'test.bam.bai'

            assert is_valid is False
            assert f'BAM index not found: {expected_bai}' in error_msg
            assert 'samtools index' in error_msg

    def test_validate_indexed_files_with_strings(self):
        """Test validation with string paths."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create test files
            fasta_path = temp_path / 'test.fasta'
            bam_path = temp_path / 'test.bam'
            fai_path = temp_path / 'test.fasta.fai'
            bai_path = temp_path / 'test.bam.bai'

            for file_path in [fasta_path, fai_path, bai_path]:
                file_path.touch()
            write_stub_bam(bam_path)

            # Pass as strings instead of Path objects
            is_valid, error_msg = validate_indexed_files(str(fasta_path), str(bam_path))

            assert is_valid is True
            assert error_msg == ''

    def test_validate_indexed_files_with_pathlib_paths(self):
        """Test validation with Path objects."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            fasta_path = temp_path / 'test.fasta'
            bam_path = temp_path / 'test.bam'
            fai_path = temp_path / 'test.fasta.fai'
            bai_path = temp_path / 'test.bam.bai'

            for file_path in [fasta_path, fai_path, bai_path]:
                file_path.touch()
            write_stub_bam(bam_path)

            is_valid, error_msg = validate_indexed_files(fasta_path, bam_path)

            assert is_valid is True
            assert error_msg == ''

    def test_sam_named_bam_is_reported_before_the_missing_index(self):
        """
        A SAM file gets the conversion recipe, not the missing-index message.

        Ordering is the whole point here. A SAM has no .bai and never will, so
        reporting the absent index first told users to run `samtools index` on
        a file that cannot be indexed.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            fasta_path = temp_path / 'test.fasta'
            fai_path = temp_path / 'test.fasta.fai'
            for file_path in [fasta_path, fai_path]:
                file_path.touch()

            # Named .bam, but the content is SAM text.
            bam_path = temp_path / 'test.bam'
            bam_path.write_bytes(b'@HD\tVN:1.6\n@SQ\tSN:contig01\tLN:560\n')

            is_valid, error_msg = validate_indexed_files(fasta_path, bam_path)

            assert is_valid is False
            assert 'samtools view -b' in error_msg
            assert 'BAM index not found' not in error_msg

    def test_validate_indexed_files_multiple_missing(self):
        """Test validation when multiple files are missing."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            fasta_path = temp_path / 'missing.fasta'
            bam_path = temp_path / 'missing.bam'

            # Don't create any files
            is_valid, error_msg = validate_indexed_files(fasta_path, bam_path)

            assert is_valid is False
            # Should fail on first missing file (FASTA)
            assert f'FASTA file not found: {fasta_path}' in error_msg

    def test_validate_indexed_files_complex_path(self):
        """Test validation with complex nested paths."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            nested_path = temp_path / 'project' / 'data' / 'alignments'
            nested_path.mkdir(parents=True)

            fasta_path = nested_path / 'genome.fasta'
            bam_path = nested_path / 'aligned_reads.bam'
            fai_path = nested_path / 'genome.fasta.fai'
            bai_path = nested_path / 'aligned_reads.bam.bai'

            for file_path in [fasta_path, fai_path, bai_path]:
                file_path.touch()
            write_stub_bam(bam_path)

            is_valid, error_msg = validate_indexed_files(fasta_path, bam_path)

            assert is_valid is True
            assert error_msg == ''

    def test_validate_indexed_files_return_types(self):
        """Test that function returns correct types."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            fasta_path = temp_path / 'test.fasta'
            bam_path = temp_path / 'missing.bam'
            fasta_path.touch()

            result = validate_indexed_files(fasta_path, bam_path)

            assert isinstance(result, tuple)
            assert len(result) == 2
            is_valid, error_msg = result
            assert isinstance(is_valid, bool)
            assert isinstance(error_msg, str)
