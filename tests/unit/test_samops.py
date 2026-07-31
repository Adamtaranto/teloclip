"""Unit tests for teloclip.samops module.

Tests SAM/BAM processing functions including CIGAR string parsing,
anchor validation, soft clip detection, and terminal position analysis.
"""

from teloclip.samops import (
    CIGARinfo,
    EnhancedStreamingSamFilter,
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
    """Test max_break and min_clip filtering in processSamlines.

    Acceptance is symmetric across both contig ends, in 1-based inclusive
    coordinates (see teloclip.overhang)::

        aln_end  = POS + reference_span - 1
        gap      = aln_start - 1        (left)  /  contig_len - aln_end  (right)
        overhang = clip_len - gap

        accepted  <=>  clip_len > 0
                       and gap <= max_break        (inclusive)
                       and overhang >= min_clip

    These tests pin the boundaries in both directions. Before the shared
    predicate existed, the left and right tests here disagreed with each other
    by a base and min_clip was not enforced on the right end at all.
    """

    def test_left_overhang_within_max_break(self):
        """Test left overhang read within max_break threshold is kept."""
        # POS 1, 50S50M: gap = 0 <= 50, overhang = 50 - 0 = 50 >= 1. Kept.
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
        # POS 51, 52S50M: gap = 50, exactly max_break, so the tolerance is
        # inclusive. overhang = 52 - 50 = 2 >= 1. Kept.
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            'read1\t0\tcontig1\t51\t30\t52S50M\t*\t0\t0\t' + 'A' * 102 + '\t*',
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
        # POS 52, 100S50M: gap = 51 > 50. Rejected on max_break, before the
        # clip is even considered.
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            'read1\t0\tcontig1\t52\t30\t100S50M\t*\t0\t0\t' + 'A' * 150 + '\t*',
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
        assert result['excluded_min_clip'] == 0

    def test_left_overhang_that_does_not_reach_contig_start(self):
        """Test left clip stopping exactly at the contig start is excluded.

        POS 51, 50S50M: gap = 50 is inside max_break, but the 50-base clip only
        reaches back to position 1 and adds nothing beyond it, so overhang = 0.
        This is a min_clip exclusion, not a max_break one.
        """
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
        assert result['excluded_max_break'] == 0
        assert result['excluded_min_clip'] == 1

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
        # POS 901, 50M52S: aln_end = 901 + 50 - 1 = 950, so gap = 1000 - 950 =
        # 50, exactly max_break. overhang = 52 - 50 = 2 >= 1. Kept.
        # This mirrors test_left_overhang_at_max_break_boundary exactly.
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            'read1\t0\tcontig1\t901\t30\t50M52S\t*\t0\t0\t' + 'A' * 102 + '\t*',
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
        # POS 899, 50M50S: aln_end = 948, gap = 1000 - 948 = 52 > 50. Excluded.
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
        # POS 901, 30M5D10M5I5M52S. The 5I consumes read bases but not
        # reference, so the reference span is 30+5+10+5 = 50 and aln_end = 950.
        # gap = 50 <= max_break, overhang = 52 - 50 = 2 >= 1. Kept.
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            'read1\t0\tcontig1\t901\t30\t30M5D10M5I5M52S\t*\t0\t0\t'
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
        # max_break=0 means the alignment must reach the terminal base exactly:
        # aln_start == 1 on the left, aln_end == contig_len on the right.
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            # POS 1: gap = 0, flush with the contig start. overhang = 2. Kept.
            'read1\t0\tcontig1\t1\t30\t2S50M\t*\t0\t0\t' + 'A' * 52 + '\t*',
            # POS 2: gap = 1 > 0. Excluded on max_break.
            'read2\t0\tcontig1\t2\t30\t4S50M\t*\t0\t0\t' + 'A' * 54 + '\t*',
            # POS 951: aln_end = 1000, gap = 0, flush with the contig end.
            # overhang = 2. Kept.
            'read3\t0\tcontig1\t951\t30\t50M2S\t*\t0\t0\t' + 'A' * 52 + '\t*',
            # POS 949: aln_end = 998, gap = 2 > 0. Excluded on max_break.
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

        assert result['keepCount'] == 2  # read1 and read3, one at each end
        assert result['excluded_max_break'] == 2  # read2 and read4

    def test_max_break_filter_with_both_overhangs(self):
        """Test max_break filter with reads having both left and right overhangs."""
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:100',
            # Read spanning the entire contig: POS 1, 27S100M27S.
            # aln_end = 100. Left: gap 0, overhang 27. Right: gap 0,
            # overhang 27. Both ends valid.
            'read1\t0\tcontig1\t1\t30\t27S100M27S\t*\t0\t0\t' + 'A' * 154 + '\t*',
            # POS 27, 30S25M30S: aln_end = 51. Left gap = 26 > 25 and right
            # gap = 49 > 25, so both ends fail on max_break.
            'read2\t0\tcontig1\t27\t30\t30S25M30S\t*\t0\t0\t' + 'A' * 85 + '\t*',
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

    def test_right_overhang_that_does_not_reach_contig_end(self):
        """Test right clip stopping exactly at the contig end is excluded.

        POS 946, 50M5S: aln_end = 995, gap = 5, and the 5-base clip only
        reaches position 1000 without passing it, so overhang = 0.

        Before the shared predicate, the right-end test was
        ``aln_end + clip >= contig_len + 1``, which reduces to overhang >= 0 and
        therefore accepted this read while ignoring min_clip entirely. The left
        end has always rejected the equivalent case.
        """
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            'read1\t0\tcontig1\t946\t30\t50M5S\t*\t0\t0\t' + 'A' * 55 + '\t*',
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
        assert result['excluded_min_clip'] == 1

    def test_right_overhang_honours_min_clip_threshold(self):
        """Test that min_clip is enforced on the right end.

        POS 946, 50M8S: aln_end = 995, gap = 5, overhang = 3. Accepted at
        min_clip 3, rejected at 4.
        """
        sam_lines = [
            '@HD\tVN:1.0\tSO:coordinate',
            '@SQ\tSN:contig1\tLN:1000',
            'read1\t0\tcontig1\t946\t30\t50M8S\t*\t0\t0\t' + 'A' * 58 + '\t*',
        ]

        from io import StringIO

        def run(min_clip):
            """
            Run processSamlines over the fixture at a given min_clip.

            Parameters
            ----------
            min_clip : int
                Minimum bases required past the contig terminus.

            Returns
            -------
            dict
                Processing counters.
            """
            return processSamlines(
                StringIO('\n'.join(sam_lines)),
                {'contig1': 1000},
                max_break=50,
                min_clip=min_clip,
                min_anchor=10,
                return_counts=True,
            )

        assert run(3)['keepCount'] == 1
        assert run(4)['keepCount'] == 0
        assert run(4)['excluded_min_clip'] == 1


class TestEnhancedStreamingSamFilterSecondary:
    """Test secondary-alignment handling in the extract streaming filter."""

    CONTIGS = {'contig1': 1000}

    @staticmethod
    def _sam_line(flag: int, name: str = 'read1') -> str:
        """
        Build a SAM line clipped at the start of contig1.

        Parameters
        ----------
        flag : int
            SAM FLAG value.
        name : str, optional
            Read name (default: 'read1').

        Returns
        -------
        str
            A tab-delimited SAM record terminated by a newline.
        """
        seq = 'A' * 150
        return (
            '\t'.join(
                [
                    name,
                    str(flag),
                    'contig1',
                    '1',
                    '60',
                    '50S100M',
                    '*',
                    '0',
                    '0',
                    seq,
                    'I' * 150,
                ]
            )
            + '\n'
        )

    def test_secondary_alignment_is_skipped_not_crashed(self):
        """Test that a secondary alignment is filtered rather than raising.

        Regression test: ``__init__`` accepted ``exclude_secondary`` but never
        stored it, while ``__iter__`` read ``self.exclude_secondary`` for every
        record carrying FLAG 256. The resulting AttributeError was not caught by
        the surrounding handler, so ``teloclip extract`` aborted on the first
        secondary alignment in the input.
        """
        lines = [self._sam_line(256, 'secondary_read')]

        results = list(
            EnhancedStreamingSamFilter(
                lines, self.CONTIGS, max_break=50, min_clip=1, min_anchor=50
            )
        )

        assert results == []

    def test_secondary_alignment_kept_when_not_excluded(self):
        """Test that secondary alignments pass through when explicitly kept."""
        lines = [self._sam_line(256, 'secondary_read')]

        results = list(
            EnhancedStreamingSamFilter(
                lines,
                self.CONTIGS,
                max_break=50,
                min_clip=1,
                min_anchor=50,
                exclude_secondary=False,
            )
        )

        assert len(results) == 1
        assert results[0]['read_name'] == 'secondary_read'
        assert results[0]['end'] == 'L'

    def test_primary_alignment_unaffected(self):
        """Test that a primary alignment is yielded regardless of the setting."""
        lines = [self._sam_line(0, 'primary_read')]

        results = list(
            EnhancedStreamingSamFilter(
                lines, self.CONTIGS, max_break=50, min_clip=1, min_anchor=50
            )
        )

        assert len(results) == 1
        assert results[0]['read_name'] == 'primary_read'


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
