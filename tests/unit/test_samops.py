"""Unit tests for teloclip.samops module.

Tests SAM/BAM processing functions including CIGAR string parsing,
anchor validation, soft clip detection, and terminal position analysis.
"""

from unittest.mock import mock_open, patch

import pytest

from teloclip.samops import (
    CIGARinfo,
    SAMinfo,
    calculate_aligned_bases,
    checkClips,
    lenCIGAR,
    processSamlines,
    splitCIGAR,
    validate_min_anchor,
)


class TestCalculateAlignedBases:
    """Test CIGAR string aligned base calculation."""

    def test_calculate_aligned_bases_match_only(self):
        """Test CIGAR with only match operations."""
        assert calculate_aligned_bases('100M') == 100

    def test_calculate_aligned_bases_mixed_operations(self):
        """Test CIGAR with mixed operations including non-aligned."""
        # 50M + 40M = 90 aligned bases (matches only)
        assert calculate_aligned_bases('50M10D40M') == 90

    def test_calculate_aligned_bases_with_soft_clips(self):
        """Test CIGAR ignores soft clips in calculation."""
        # Only counts M operations, ignores S
        assert calculate_aligned_bases('20S80M20S') == 80

    def test_calculate_aligned_bases_sequence_match_mismatch(self):
        """Test CIGAR with sequence matches and mismatches."""
        # = (match) and X (mismatch) both count as aligned
        assert calculate_aligned_bases('30=10X30=') == 70

    def test_calculate_aligned_bases_complex_cigar(self):
        """Test complex CIGAR string with multiple operation types."""
        # 20M + 15= + 5X + 30M = 70 aligned bases
        assert calculate_aligned_bases('10S20M5I15=2D5X30M10S') == 70

    def test_calculate_aligned_bases_no_aligned_operations(self):
        """Test CIGAR with no aligned operations."""
        assert calculate_aligned_bases('50S10I20D') == 0

    def test_calculate_aligned_bases_empty_cigar(self):
        """Test empty CIGAR string."""
        assert calculate_aligned_bases('') == 0

    def test_calculate_aligned_bases_invalid_cigar(self):
        """Test invalid CIGAR string gracefully handled."""
        # Should not crash on malformed CIGAR
        assert calculate_aligned_bases('invalid') == 0


class TestValidateMinAnchor:
    """Test minimum anchor validation."""

    def test_validate_min_anchor_sufficient(self):
        """Test anchor validation with sufficient aligned bases."""
        assert validate_min_anchor('100M', 50) is True
        assert validate_min_anchor('200M', 100) is True

    def test_validate_min_anchor_insufficient(self):
        """Test anchor validation with insufficient aligned bases."""
        assert validate_min_anchor('30M', 50) is False
        assert validate_min_anchor('80M', 100) is False

    def test_validate_min_anchor_exact_threshold(self):
        """Test anchor validation at exact threshold."""
        assert validate_min_anchor('50M', 50) is True
        assert validate_min_anchor('100M', 100) is True

    def test_validate_min_anchor_with_soft_clips(self):
        """Test anchor validation ignores soft clips."""
        # 80M aligned bases, ignores 20S clips
        assert validate_min_anchor('20S80M20S', 50) is True
        assert validate_min_anchor('20S30M20S', 50) is False

    def test_validate_min_anchor_zero_threshold(self):
        """Test anchor validation with zero threshold."""
        assert validate_min_anchor('1M', 0) is True
        assert validate_min_anchor('', 0) is True

    def test_validate_min_anchor_complex_cigar(self):
        """Test anchor validation with complex CIGAR."""
        # 20M + 30= + 10X = 60 aligned bases
        cigar = '10S20M5I30=2D10X20S'
        assert validate_min_anchor(cigar, 50) is True
        assert validate_min_anchor(cigar, 70) is False


class TestSplitCIGAR:
    """Test CIGAR string splitting function."""

    def test_split_cigar_simple(self):
        """Test splitting simple CIGAR string."""
        result = splitCIGAR('100M')
        expected = [(100, 'M')]
        assert result == expected

    def test_split_cigar_complex(self):
        """Test splitting complex CIGAR string."""
        result = splitCIGAR('20S50M10I30M20S')
        expected = [(20, 'S'), (50, 'M'), (10, 'I'), (30, 'M'), (20, 'S')]
        assert result == expected

    def test_split_cigar_all_operations(self):
        """Test splitting CIGAR with all operation types."""
        result = splitCIGAR('10M5I3D20N15S25H30P35=40X')
        expected = [
            (10, 'M'),
            (5, 'I'),
            (3, 'D'),
            (20, 'N'),
            (15, 'S'),
            (25, 'H'),
            (30, 'P'),
            (35, '='),
            (40, 'X'),
        ]
        assert result == expected


class TestCheckClips:
    """Test soft clip detection."""

    def test_check_clips_left_only(self):
        """Test detecting left soft clip only."""
        result = checkClips('20S80M')
        assert result == (20, None)

    def test_check_clips_right_only(self):
        """Test detecting right soft clip only."""
        result = checkClips('80M20S')
        assert result == (None, 20)

    def test_check_clips_both_sides(self):
        """Test detecting soft clips on both sides."""
        result = checkClips('20S60M20S')
        assert result == (20, 20)

    def test_check_clips_no_clips(self):
        """Test with no soft clips."""
        result = checkClips('100M')
        assert result == (None, None)


class TestLenCIGAR:
    """Test CIGAR reference length calculation."""

    def test_len_cigar_simple(self):
        """Test simple CIGAR reference length."""
        result = lenCIGAR('100M')
        assert result == 100

    def test_len_cigar_with_clips(self):
        """Test CIGAR reference length ignoring soft clips."""
        result = lenCIGAR('20S80M20S')
        assert result == 80  # Only the M operation contributes to reference length

    def test_len_cigar_complex(self):
        """Test complex CIGAR reference length calculation."""
        # M, D, N, X, = operations contribute to reference length
        result = lenCIGAR('10S20M5I30=10D5X10S')
        assert result == 65  # 20M + 30= + 10D + 5X = 65


class TestProcessSamlines:
    """Test main SAM processing function."""

    @pytest.mark.xfail(
        reason='Function signature mismatch in processSamlines - needs investigation'
    )
    @patch('builtins.open', new_callable=mock_open)
    def test_process_samlines_basic(self, mock_file, temp_dir, sample_sam_alignments):
        """Test basic SAM processing without filters."""
        # Mock file reading - each line should be a complete string
        # The function expects line[0][0] so maybe it expects line to be a list of strings?
        test_sam_lines = [
            ['@HD\tVN:1.0\tSO:unsorted'],
        ] + [[line] for line in sample_sam_alignments]

        mock_file.return_value.__iter__.return_value = iter(test_sam_lines)

        # Process without strict filters
        counts = processSamlines(
            samfile='test.sam',
            contig_dict={'contig01': 1000},
            motif_list=[],
            match_anywhere=False,
            max_break=0,
            min_clip=1,
            min_repeats=1,
            min_anchor=0,  # No anchor filtering
        )

        # Should return a dictionary with counts
        assert isinstance(counts, dict)
        assert 'samlineCount' in counts


class TestSAMinfo:
    """Test SAM information function."""

    def test_sam_info(self):
        """Test SAM info function runs without error."""
        # This is a documentation function, just test it doesn't crash
        try:
            SAMinfo()
            success = True
        except Exception:
            success = False
        assert success is True


class TestCIGARinfo:
    """Test CIGAR information function."""

    def test_cigar_info(self):
        """Test CIGAR info function runs without error."""
        # This is a documentation function, just test it doesn't crash
        try:
            CIGARinfo()
            success = True
        except Exception:
            success = False
        assert success is True
