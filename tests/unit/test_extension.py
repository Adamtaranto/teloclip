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

    def test_validate_accepts_trimmed_extension(self):
        """Test validation of a net-positive extension that trimmed bases."""
        original = 'ATCGATCGATCG'
        overhang = OverhangInfo('GGGG', 4, 3, 100, 'read1', True, 4, 96, 'test_contig')
        # First 2 bases trimmed, 4-base overhang prepended: net +2.
        extended = 'GGGGCGATCGATCG'

        assert validate_extension(original, extended, overhang, trim_length=2) is True

    def test_validate_rejects_net_shortening_after_trim(self):
        """Test validation of an extension that trims more than it adds.

        The length comparison must be made against the untrimmed original.
        Comparing against the trimmed remainder makes this trivially true, which
        is how net-shortening extensions previously went undetected.
        """
        original = 'ATCGATCGATCG'
        overhang = OverhangInfo('GGG', 3, 10, 12, 'read1', True, 3, 96, 'test_contig')
        # 9 bases trimmed, 3-base overhang prepended: net -6.
        extended = 'GGG' + original[9:]

        assert validate_extension(original, extended, overhang, trim_length=9) is False


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

    def test_apply_extension_zero_gain_rejected(self):
        """Test that an overhang adding nothing is rejected."""
        contig_seq = 'ATCGATCGATCG'
        # An empty overhang at a flush alignment: no trim, nothing added.
        overhang = OverhangInfo('', 0, 1, 100, 'read1', True, 0, 100, 'test_contig')

        with pytest.raises(ValueError, match=r'would change its length by \+0 bp'):
            apply_contig_extension(contig_seq, overhang, len(contig_seq))

    def test_apply_left_extension_that_would_shorten_is_rejected(self):
        """Test that a short left clip anchored inside the contig is rejected.

        Regression test for the core defect this guard exists to prevent. A read
        aligning from position 10 of a 12 bp contig with only a 3 bp clip causes
        9 bases to be trimmed and 3 to be added, shortening the contig by 6.
        Before the guard, ``validate_extension`` compared the result against the
        already-trimmed sequence, so this passed silently.
        """
        contig_seq = 'ATCGATCGATCG'
        overhang = OverhangInfo('GGG', 3, 10, 12, 'read1', True, 3, 96, 'test_contig')

        with pytest.raises(ValueError, match=r'would change its length by -6 bp'):
            apply_contig_extension(contig_seq, overhang, len(contig_seq))

    def test_apply_right_extension_that_would_shorten_is_rejected(self):
        """Test that a short right clip anchored inside the contig is rejected."""
        contig_seq = 'ATCGATCGATCG'
        # Alignment ends at position 3 of a 12 bp contig: 9 bases trimmed.
        overhang = OverhangInfo('GGG', 3, 1, 3, 'read1', False, 3, 96, 'test_contig')

        with pytest.raises(ValueError, match=r'would change its length by -6 bp'):
            apply_contig_extension(contig_seq, overhang, len(contig_seq))

    def test_apply_extension_with_zero_net_gain_rejected(self):
        """Test that a clip exactly reaching the terminus is rejected.

        Trimming 4 bases and adding 4 leaves the contig the same length, which
        is not an extension.
        """
        contig_seq = 'ATCGATCGATCG'
        overhang = OverhangInfo('GGGG', 4, 5, 12, 'read1', True, 4, 96, 'test_contig')

        with pytest.raises(ValueError, match=r'would change its length by \+0 bp'):
            apply_contig_extension(contig_seq, overhang, len(contig_seq))

    def test_min_net_gain_is_configurable(self):
        """Test that the required net gain can be raised above the default."""
        contig_seq = 'ATCGATCGATCG'
        # Trims 2, adds 4: a net gain of exactly 2.
        overhang = OverhangInfo('GGGG', 4, 3, 100, 'read1', True, 4, 96, 'test_contig')

        _, ext_info = apply_contig_extension(
            contig_seq, overhang, len(contig_seq), min_net_gain=2
        )
        assert ext_info['net_gain'] == 2

        with pytest.raises(ValueError, match='at least 3 bp is required'):
            apply_contig_extension(
                contig_seq, overhang, len(contig_seq), min_net_gain=3
            )

    def test_net_gain_reconciles_with_lengths(self):
        """Test that original + net_gain == final for a trimmed extension."""
        contig_seq = 'ATCGATCGATCG'
        overhang = OverhangInfo('GGGG', 4, 3, 100, 'read1', True, 4, 96, 'test_contig')

        _, ext_info = apply_contig_extension(contig_seq, overhang, len(contig_seq))

        assert ext_info['net_gain'] == 2
        assert (
            ext_info['original_length'] + ext_info['net_gain']
            == ext_info['final_length']
        )
        # And the same identity expressed through the raw clip and trim.
        assert (
            ext_info['original_length']
            + ext_info['overhang_length']
            - ext_info['trim_length']
            == ext_info['final_length']
        )

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
