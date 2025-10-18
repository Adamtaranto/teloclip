"""
Tests for anchor validation functions in samops module.
"""

from teloclip.samops import (
    calculate_aligned_bases,
    checkClips,
    lenCIGAR,
    splitCIGAR,
    validate_min_anchor,
)


class TestCalculateAlignedBases:
    """Test calculate_aligned_bases function."""

    def test_simple_match(self):
        """Test simple match operations."""
        assert calculate_aligned_bases('100M') == 100
        assert calculate_aligned_bases('50M') == 50

    def test_sequence_match_and_mismatch(self):
        """Test sequence match and mismatch operations."""
        assert calculate_aligned_bases('50=') == 50  # sequence match
        assert calculate_aligned_bases('25X') == 25  # mismatch
        assert calculate_aligned_bases('30=20X') == 50  # both

    def test_mixed_match_operations(self):
        """Test combinations of M, =, and X operations."""
        assert calculate_aligned_bases('30M20=10X') == 60
        assert calculate_aligned_bases('10M5X15=') == 30

    def test_with_soft_clips(self):
        """Test that soft clips are excluded from aligned bases count."""
        assert calculate_aligned_bases('20S100M') == 100
        assert calculate_aligned_bases('100M30S') == 100
        assert calculate_aligned_bases('10S50M20S') == 50

    def test_with_insertions_deletions(self):
        """Test that insertions and deletions are excluded."""
        assert calculate_aligned_bases('50M10I40M') == 90  # I excluded
        assert calculate_aligned_bases('50M5D40M') == 90  # D excluded
        assert calculate_aligned_bases('30M10I5D20M') == 50  # both excluded

    def test_with_hard_clips(self):
        """Test that hard clips are excluded."""
        assert calculate_aligned_bases('15H100M') == 100
        assert calculate_aligned_bases('100M25H') == 100
        assert calculate_aligned_bases('10H50M15H') == 50

    def test_with_padding_and_splicing(self):
        """Test that padding and splicing are excluded."""
        assert calculate_aligned_bases('50M10P40M') == 90  # P excluded
        assert calculate_aligned_bases('50M20N40M') == 90  # N excluded

    def test_complex_cigar(self):
        """Test complex CIGAR strings with multiple operations."""
        # 20S30M10I5D25M15S: only 30M + 25M = 55 aligned bases
        assert calculate_aligned_bases('20S30M10I5D25M15S') == 55

        # 10H5S40M20I10D30=20X5S10H: only 40M + 30= + 20X = 90 aligned bases
        assert calculate_aligned_bases('10H5S40M20I10D30=20X5S10H') == 90

    def test_only_clips_and_indels(self):
        """Test CIGAR strings with no aligned bases."""
        assert calculate_aligned_bases('50S') == 0
        assert calculate_aligned_bases('20H30S') == 0
        assert calculate_aligned_bases('10S20I30S') == 0
        assert calculate_aligned_bases('100I') == 0


class TestValidateMinAnchor:
    """Test validate_min_anchor function."""

    def test_meets_requirement(self):
        """Test alignments that meet min_anchor requirement."""
        assert validate_min_anchor('100M', 50)
        assert validate_min_anchor('100M', 100)
        assert validate_min_anchor('200M', 150)

    def test_fails_requirement(self):
        """Test alignments that fail min_anchor requirement."""
        assert not validate_min_anchor('50M', 100)
        assert not validate_min_anchor('25M', 50)
        assert not validate_min_anchor('99M', 100)

    def test_exact_requirement(self):
        """Test alignments that exactly meet requirement."""
        assert validate_min_anchor('100M', 100)
        assert validate_min_anchor('50M', 50)
        assert validate_min_anchor('1M', 1)

    def test_with_soft_clips(self):
        """Test that soft clips don't count toward anchor."""
        assert validate_min_anchor('100S50M', 50)  # meets exactly
        assert not validate_min_anchor('100S50M', 51)  # fails by 1
        assert validate_min_anchor('50M100S', 50)  # meets exactly
        assert validate_min_anchor('20S30M40S', 30)  # meets exactly

    def test_with_indels(self):
        """Test with insertions and deletions."""
        assert validate_min_anchor('30M10I20M', 50)  # 30+20=50, meets exactly
        assert not validate_min_anchor('30M10I20M', 51)  # fails by 1
        assert validate_min_anchor('30M5D20M', 50)  # 30+20=50, meets exactly

    def test_complex_scenarios(self):
        """Test complex real-world scenarios."""
        # Read with terminal soft clip but sufficient anchor
        assert validate_min_anchor('600M50S', 500)
        assert not validate_min_anchor('600M50S', 650)  # clips don't count

        # Read with leading soft clip
        assert validate_min_anchor('50S600M', 500)
        assert not validate_min_anchor('50S400M', 500)

        # Read with clips on both ends
        assert validate_min_anchor('50S500M100S', 400)
        assert not validate_min_anchor('50S500M100S', 600)  # clips excluded

    def test_zero_anchor_requirement(self):
        """Test with zero anchor requirement."""
        assert validate_min_anchor('100S', 0)
        assert validate_min_anchor('1M', 0)
        assert validate_min_anchor('50S50M', 0)


class TestAnchorIntegrationWithExistingFunctions:
    """Test anchor functions work correctly with existing CIGAR parsing."""

    def test_consistency_with_split_cigar(self):
        """Test that our anchor calculation is consistent with splitCIGAR."""
        test_cigars = [
            '100M',
            '50M50S',
            '20S80M',
            '30M10I20D40M',
            '10H20S50M30=20X40S15H',
        ]

        for cigar in test_cigars:
            # Our function should give same result as manual calculation
            cigar_list = splitCIGAR(cigar)
            manual_count = sum(
                length for length, op in cigar_list if op in {'M', '=', 'X'}
            )
            assert calculate_aligned_bases(cigar) == manual_count

    def test_anchor_vs_alignment_length(self):
        """Test difference between anchor length and total alignment length."""
        # CIGAR with deletions - alignment length includes D, anchor length doesn't
        cigar = '50M10D40M'
        anchor_length = calculate_aligned_bases(cigar)  # 50 + 40 = 90
        alignment_length = lenCIGAR(cigar)  # 50 + 10 + 40 = 100

        assert anchor_length == 90
        assert alignment_length == 100
        assert anchor_length < alignment_length

        # CIGAR with only matches - should be equal
        cigar = '100M'
        anchor_length = calculate_aligned_bases(cigar)
        alignment_length = lenCIGAR(cigar)

        assert anchor_length == 100
        assert alignment_length == 100
        assert anchor_length == alignment_length

    def test_clips_detection_consistency(self):
        """Test that clips are detected consistently."""
        test_cases = [
            ('50S100M', True, False),  # left clip only
            ('100M50S', False, True),  # right clip only
            ('20S80M30S', True, True),  # both clips
            ('100M', False, False),  # no clips
        ]

        for cigar, expect_left, expect_right in test_cases:
            left_clip, right_clip = checkClips(cigar)
            has_left = left_clip is not None
            has_right = right_clip is not None

            assert has_left == expect_left, f'Left clip detection failed for {cigar}'
            assert has_right == expect_right, f'Right clip detection failed for {cigar}'

            # If there are clips, anchor should be less than total read length
            if has_left or has_right:
                anchor_bases = calculate_aligned_bases(cigar)
                # We can't easily calculate read length without more parsing,
                # but we can verify anchor excludes clip operations
                cigar_list = splitCIGAR(cigar)
                total_ops = sum(length for length, op in cigar_list)
                clip_bases = sum(length for length, op in cigar_list if op == 'S')
                if clip_bases > 0:
                    assert anchor_bases <= total_ops - clip_bases
