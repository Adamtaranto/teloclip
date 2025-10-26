"""Unit tests for teloclip.samops module.

Tests SAM/BAM processing functions including CIGAR string parsing,
anchor validation, soft clip detection, and terminal position analysis.
"""

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

    def test_process_samlines_basic(self, sample_sam_alignments):
        """Test basic SAM processing without filters."""
        # Create test SAM lines as strings (not nested lists)
        test_sam_lines = ['@HD\tVN:1.0\tSO:unsorted'] + sample_sam_alignments

        # Process without strict filters, requesting counts for testing
        counts = processSamlines(
            samfile=test_sam_lines,
            contig_dict={
                'contig01': 1000,
                'contig02': 2000,
                'contig03': 1500,
            },
            motif_list=[],
            match_anywhere=False,
            max_break=0,
            min_clip=1,
            min_repeats=1,
            min_anchor=0,  # No anchor filtering
            return_counts=True,  # Request counts for testing
        )

        # Should return a dictionary with counts
        assert isinstance(counts, dict)
        assert 'samlineCount' in counts
        assert 'keepCount' in counts
        assert 'motifCount' in counts
        assert 'removeCount' in counts
        assert 'anchorFilteredCount' in counts
        assert 'bothCount' in counts

        # Should have processed some SAM records (not counting header)
        assert counts['samlineCount'] > 0


class TestProcessSamlinesMaxBreakFilter:
    """Test max_break filter application in processSamlines.

    These tests verify that the max_break filter is correctly applied to both
    left and right overhangs in SAM alignments:

    Left overhangs:
    - Position (SAM_POS) must be <= max_break from contig start (1-based)
    - Soft clip length must be >= (position + min_clip) to extend past contig start

    Right overhangs:
    - Alignment end must be <= max_break from contig end
    - Soft clip must extend >= 1 base past contig end (alnEnd + clipLen >= contigLen + 1)

    The tests verify boundary conditions, exclusions, and complex CIGAR strings.
    """

    def test_left_overhang_within_max_break(self):
        """Test left overhang read within max_break threshold is kept."""
        # Read at position 1 with 50bp soft clip, max_break=50
        # Should be kept since position 1 <= max_break (50)
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            'read1\t0\tcontig1\t1\t30\t50S50M\t*\t0\t0\t' + 'A' * 100 + '\t*',
        ]

        from io import StringIO

        samfile = StringIO('\n'.join(sam_lines))
        contig_dict = {'contig1': 1000}

        result = processSamlines(
            samfile,
            contig_dict,
            max_break=50,
            min_clip=1,
            min_anchor=10,
            return_counts=True,
        )

        assert result['keepCount'] == 1
        assert result['excluded_max_break'] == 0

    def test_left_overhang_at_max_break_boundary(self):
        """Test left overhang read exactly at max_break threshold is kept."""
        # Read at position 50 with 52bp soft clip, max_break=50, min_clip=1
        # Conditions: pos <= max_break (50 <= 50) AND leftClipLen >= (pos + min_clip) (52 >= 51)
        # Should be kept since both conditions are met
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            'read1\t0\tcontig1\t50\t30\t52S50M\t*\t0\t0\t' + 'A' * 102 + '\t*',
        ]

        from io import StringIO

        samfile = StringIO('\n'.join(sam_lines))
        contig_dict = {'contig1': 1000}

        result = processSamlines(
            samfile,
            contig_dict,
            max_break=50,
            min_clip=1,
            min_anchor=10,
            return_counts=True,
        )

        assert result['keepCount'] == 1
        assert result['excluded_max_break'] == 0

    def test_left_overhang_exceeds_max_break(self):
        """Test left overhang read beyond max_break threshold is excluded."""
        # Read at position 51 with 50bp soft clip, max_break=50
        # Should be excluded since position 51 > max_break (50)
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            'read1\t0\tcontig1\t51\t30\t50S50M\t*\t0\t0\t' + 'A' * 100 + '\t*',
        ]

        from io import StringIO

        samfile = StringIO('\n'.join(sam_lines))
        contig_dict = {'contig1': 1000}

        result = processSamlines(
            samfile,
            contig_dict,
            max_break=50,
            min_clip=1,
            min_anchor=10,
            return_counts=True,
        )

        assert result['keepCount'] == 0
        assert result['excluded_max_break'] == 1

    def test_right_overhang_within_max_break(self):
        """Test right overhang read within max_break threshold is kept."""
        # Contig length 1000, read alignment ends at position 950 + 50 = 1000
        # Distance from contig end: 1000 - 1000 = 0 <= max_break (50)
        # Should be kept
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            'read1\t0\tcontig1\t950\t30\t50M50S\t*\t0\t0\t' + 'A' * 100 + '\t*',
        ]

        from io import StringIO

        samfile = StringIO('\n'.join(sam_lines))
        contig_dict = {'contig1': 1000}

        result = processSamlines(
            samfile,
            contig_dict,
            max_break=50,
            min_clip=1,
            min_anchor=10,
            return_counts=True,
        )

        assert result['keepCount'] == 1
        assert result['excluded_max_break'] == 0

    def test_right_overhang_at_max_break_boundary(self):
        """Test right overhang read exactly at max_break threshold is kept."""
        # Contig length 1000, read alignment ends at position 900 + 50 = 950
        # Distance from contig end: 1000 - 950 = 50 <= max_break (50)
        # Overhang condition: alnEnd + rightClipLen >= ContigLen + 1 → 950 + 52 >= 1001 → 1002 >= 1001 ✓
        # Should be kept
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            'read1\t0\tcontig1\t900\t30\t50M52S\t*\t0\t0\t' + 'A' * 102 + '\t*',
        ]

        from io import StringIO

        samfile = StringIO('\n'.join(sam_lines))
        contig_dict = {'contig1': 1000}

        result = processSamlines(
            samfile,
            contig_dict,
            max_break=50,
            min_clip=1,
            min_anchor=10,
            return_counts=True,
        )

        assert result['keepCount'] == 1
        assert result['excluded_max_break'] == 0

    def test_right_overhang_exceeds_max_break(self):
        """Test right overhang read beyond max_break threshold is excluded."""
        # Contig length 1000, read alignment ends at position 899 + 50 = 949
        # Distance from contig end: 1000 - 949 = 51 > max_break (50)
        # Should be excluded
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            'read1\t0\tcontig1\t899\t30\t50M50S\t*\t0\t0\t' + 'A' * 100 + '\t*',
        ]

        from io import StringIO

        samfile = StringIO('\n'.join(sam_lines))
        contig_dict = {'contig1': 1000}

        result = processSamlines(
            samfile,
            contig_dict,
            max_break=50,
            min_clip=1,
            min_anchor=10,
            return_counts=True,
        )

        assert result['keepCount'] == 0
        assert result['excluded_max_break'] == 1

    def test_right_overhang_with_gaps_in_cigar(self):
        """Test right overhang with deletions/insertions in CIGAR string."""
        # Contig length 1000, read starts at 900
        # CIGAR: 30M5D10M5I5M = alignment length on reference = 30+5+10+5 = 50
        # Alignment end: 900 + 50 = 950
        # Distance from contig end: 1000 - 950 = 50 <= max_break (50)
        # Overhang condition: alnEnd + rightClipLen >= ContigLen + 1 → 950 + 52 >= 1001 → 1002 >= 1001 ✓
        # Should be kept
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            'read1\t0\tcontig1\t900\t30\t30M5D10M5I5M52S\t*\t0\t0\t'
            + 'A' * 102
            + '\t*',
        ]

        from io import StringIO

        samfile = StringIO('\n'.join(sam_lines))
        contig_dict = {'contig1': 1000}

        result = processSamlines(
            samfile,
            contig_dict,
            max_break=50,
            min_clip=1,
            min_anchor=10,
            return_counts=True,
        )

        assert result['keepCount'] == 1
        assert result['excluded_max_break'] == 0

    def test_zero_max_break_filter(self):
        """Test max_break=0 only allows reads at exact contig boundaries."""
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            # Left overhang at position 1 (not allowed with max_break=0 since 1 > 0)
            'read1\t0\tcontig1\t1\t30\t2S50M\t*\t0\t0\t' + 'A' * 52 + '\t*',
            # Left overhang at position 2 (not allowed with max_break=0)
            'read2\t0\tcontig1\t2\t30\t4S50M\t*\t0\t0\t' + 'A' * 54 + '\t*',
            # Right overhang ending exactly at contig end: 950+50=1000, clip 2 → 1000+2=1002 >= 1001 ✓
            'read3\t0\tcontig1\t950\t30\t50M2S\t*\t0\t0\t' + 'A' * 52 + '\t*',
            # Right overhang ending 1bp before contig end: 949+50=999, distance=1000-999=1 > max_break(0)
            'read4\t0\tcontig1\t949\t30\t50M52S\t*\t0\t0\t' + 'A' * 102 + '\t*',
        ]

        from io import StringIO

        samfile = StringIO('\n'.join(sam_lines))
        contig_dict = {'contig1': 1000}

        result = processSamlines(
            samfile,
            contig_dict,
            max_break=0,
            min_clip=1,
            min_anchor=10,
            return_counts=True,
        )

        assert result['keepCount'] == 1  # Only read3
        assert result['excluded_max_break'] == 3  # read1, read2, and read4

    def test_max_break_filter_with_both_overhangs(self):
        """Test max_break filter with reads having both left and right overhangs."""
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:100',
            # Read spanning entire contig: pos 1, length 100, both overhangs valid
            # Left: pos 1 <= 25, clip 27 >= (1+1)=2 ✓
            # Right: alnEnd=101, (100-101)=-1 <= 25 ✓, 101+27=128 >= 101 ✓
            'read1\t0\tcontig1\t1\t30\t27S100M27S\t*\t0\t0\t' + 'A' * 154 + '\t*',
            # Read with BOTH overhangs exceeding max_break to ensure exclusion
            # Left: pos 50 > 25, Right: alnEnd=100, (100-100)=0 <= 25 but pos 50 > 25 for left
            # To make right also fail: use position 26, alnEnd = 26+25 = 51, (100-51)=49 > 25
            'read2\t0\tcontig1\t26\t30\t30S25M30S\t*\t0\t0\t' + 'A' * 85 + '\t*',
        ]

        from io import StringIO

        samfile = StringIO('\n'.join(sam_lines))
        contig_dict = {'contig1': 100}

        result = processSamlines(
            samfile,
            contig_dict,
            max_break=25,
            min_clip=1,
            min_anchor=10,
            return_counts=True,
        )

        assert result['keepCount'] == 1  # Only read1 should be kept
        assert result['bothCount'] == 1  # read1 spans entire contig
        assert result['excluded_max_break'] == 1  # read2 excluded


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
