"""
Unit tests for teloclip.extension module.

Tests for contig extension algorithms and position calculation.
"""

import pytest

from teloclip.analysis import OverhangInfo
from teloclip.extension import (
    apply_contig_extension,
    calculate_extension_position,
    extend_contig,
    trim_contig_end,
    validate_extension,
)


class TestCalculateExtensionPosition:
    """Test the calculate_extension_position function."""

    def test_left_extension_exact_start(self):
        """Test left extension when alignment starts exactly at position 1."""
        position, trim_length = calculate_extension_position(
            alignment_pos=1, alignment_end=100, contig_length=1000, is_left=True
        )

        assert position == 0
        assert trim_length == 0

    def test_left_extension_with_gap(self):
        """Test left extension when alignment starts after position 1."""
        position, trim_length = calculate_extension_position(
            alignment_pos=5, alignment_end=100, contig_length=1000, is_left=True
        )

        assert position == 0
        assert trim_length == 4  # Need to trim 4 bases (positions 1-4)

    def test_right_extension_exact_end(self):
        """Test right extension when alignment ends exactly at contig end."""
        position, trim_length = calculate_extension_position(
            alignment_pos=900, alignment_end=1000, contig_length=1000, is_left=False
        )

        assert position == 1000
        assert trim_length == 0

    def test_right_extension_with_gap(self):
        """Test right extension when alignment ends before contig end."""
        position, trim_length = calculate_extension_position(
            alignment_pos=900, alignment_end=995, contig_length=1000, is_left=False
        )

        assert position == 995
        assert trim_length == 5  # Need to trim 5 bases from end


class TestTrimContigEnd:
    """Test the trim_contig_end function."""

    def test_trim_left_end(self):
        """Test trimming from left end of contig."""
        sequence = 'ATCGATCGATCG'
        result = trim_contig_end(sequence, 4, is_left_end=True)

        assert result == 'ATCGATCG'

    def test_trim_right_end(self):
        """Test trimming from right end of contig."""
        sequence = 'ATCGATCGATCG'
        result = trim_contig_end(sequence, 4, is_left_end=False)

        assert result == 'ATCGATCG'

    def test_trim_zero_length(self):
        """Test trimming zero bases."""
        sequence = 'ATCGATCGATCG'
        result = trim_contig_end(sequence, 0, is_left_end=True)

        assert result == sequence

    def test_trim_negative_length(self):
        """Test trimming negative length (should do nothing)."""
        sequence = 'ATCGATCGATCG'
        result = trim_contig_end(sequence, -5, is_left_end=True)

        assert result == sequence


class TestExtendContig:
    """Test the extend_contig function."""

    def test_extend_left_end(self):
        """Test extending left end of contig."""
        sequence = 'ATCGATCGATCG'
        overhang = OverhangInfo('GGGG', 4, 1, 100, 'read1', True, 4, 96, 'test_contig')

        result = extend_contig(sequence, overhang, 0, is_left_end=True)

        assert result == 'GGGGATCGATCGATCG'

    def test_extend_right_end(self):
        """Test extending right end of contig."""
        sequence = 'ATCGATCGATCG'
        overhang = OverhangInfo(
            'AAAA', 4, 995, 1000, 'read1', False, 4, 96, 'test_contig'
        )

        result = extend_contig(sequence, overhang, len(sequence), is_left_end=False)

        assert result == 'ATCGATCGATCGAAAA'


class TestValidateExtension:
    """Test the validate_extension function."""

    def test_validate_left_extension_success(self):
        """Test validation of successful left extension."""
        original = 'ATCGATCGATCG'
        overhang = OverhangInfo('GGGG', 4, 1, 100, 'read1', True, 4, 96, 'test_contig')
        extended = 'GGGGATCGATCGATCG'

        result = validate_extension(original, extended, overhang)

        assert result is True

    def test_validate_right_extension_success(self):
        """Test validation of successful right extension."""
        original = 'ATCGATCGATCG'
        overhang = OverhangInfo(
            'AAAA', 4, 995, 1000, 'read1', False, 4, 96, 'test_contig'
        )
        extended = 'ATCGATCGATCGAAAA'

        result = validate_extension(original, extended, overhang)

        assert result is True

    def test_validate_left_extension_failure(self):
        """Test validation of failed left extension."""
        original = 'ATCGATCGATCG'
        overhang = OverhangInfo('GGGG', 4, 1, 100, 'read1', True, 4, 96, 'test_contig')
        extended = 'ATCGATCGATCGAAAA'  # Wrong extension

        result = validate_extension(original, extended, overhang)

        assert result is False

    def test_validate_right_extension_failure(self):
        """Test validation of failed right extension."""
        original = 'ATCGATCGATCG'
        overhang = OverhangInfo(
            'AAAA', 4, 995, 1000, 'read1', False, 4, 96, 'test_contig'
        )
        extended = 'GGGGATCGATCGATCG'  # Wrong extension

        result = validate_extension(original, extended, overhang)

        assert result is False

    def test_validate_no_length_increase(self):
        """Test validation when extended sequence is not longer."""
        original = 'ATCGATCGATCG'
        overhang = OverhangInfo(
            'AAAA', 4, 995, 1000, 'read1', False, 4, 96, 'test_contig'
        )
        extended = 'ATCGATCG'  # Shorter than original

        result = validate_extension(original, extended, overhang)

        assert result is False


class TestApplyContigExtension:
    """Test the apply_contig_extension function."""

    def test_apply_left_extension_no_trim(self):
        """Test applying left extension without trimming."""
        contig_seq = 'ATCGATCGATCG'
        overhang = OverhangInfo('GGGG', 4, 1, 100, 'read1', True, 4, 96, 'test_contig')

        extended_seq, ext_info = apply_contig_extension(
            contig_seq, overhang, len(contig_seq)
        )

        assert extended_seq == 'GGGGATCGATCGATCG'
        assert ext_info['overhang_length'] == 4
        assert ext_info['trim_length'] == 0
        assert ext_info['original_length'] == 12
        assert ext_info['final_length'] == 16
        assert ext_info['is_left'] is True

    def test_apply_right_extension_no_trim(self):
        """Test applying right extension without trimming."""
        contig_seq = 'ATCGATCGATCG'
        overhang = OverhangInfo('AAAA', 4, 12, 12, 'read1', False, 4, 96, 'test_contig')

        extended_seq, ext_info = apply_contig_extension(
            contig_seq, overhang, len(contig_seq)
        )

        assert extended_seq == 'ATCGATCGATCGAAAA'
        assert ext_info['overhang_length'] == 4
        assert ext_info['trim_length'] == 0
        assert ext_info['is_left'] is False

    def test_apply_left_extension_with_trim(self):
        """Test applying left extension with trimming."""
        contig_seq = 'ATCGATCGATCG'
        overhang = OverhangInfo('GGGG', 4, 3, 100, 'read1', True, 4, 96, 'test_contig')

        extended_seq, ext_info = apply_contig_extension(
            contig_seq, overhang, len(contig_seq)
        )

        # Should trim first 2 bases and add overhang
        assert extended_seq == 'GGGGCGATCGATCG'
        assert ext_info['trim_length'] == 2

    def test_apply_right_extension_with_trim(self):
        """Test applying right extension with trimming."""
        contig_seq = 'ATCGATCGATCG'
        overhang = OverhangInfo(
            'AAAA', 4, 900, 10, 'read1', False, 4, 96, 'test_contig'
        )

        extended_seq, ext_info = apply_contig_extension(
            contig_seq, overhang, len(contig_seq)
        )

        # Should trim last 2 bases and add overhang
        assert extended_seq == 'ATCGATCGATAAAA'
        assert ext_info['trim_length'] == 2

    def test_apply_extension_validation_failure(self):
        """Test that validation failure raises exception."""
        contig_seq = 'ATCGATCGATCG'
        # Create an overhang that would cause validation to fail
        # This is a bit artificial but tests the error path
        overhang = OverhangInfo('', 0, 1, 100, 'read1', True, 0, 100, 'test_contig')

        with pytest.raises(ValueError, match='Extension validation failed'):
            apply_contig_extension(contig_seq, overhang, len(contig_seq))

    def test_extension_info_completeness(self):
        """Test that extension info contains all expected fields."""
        contig_seq = 'ATCGATCGATCG'
        overhang = OverhangInfo('GGGG', 4, 1, 100, 'read1', True, 4, 96, 'test_contig')

        extended_seq, ext_info = apply_contig_extension(
            contig_seq, overhang, len(contig_seq)
        )

        required_fields = [
            'overhang_length',
            'trim_length',
            'extension_position',
            'original_length',
            'final_length',
            'read_name',
            'is_left',
        ]

        for field in required_fields:
            assert field in ext_info

        assert ext_info['read_name'] == 'read1'
