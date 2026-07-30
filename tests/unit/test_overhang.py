"""
Unit tests for the canonical overhang module.

These tests pin the coordinate convention and the acceptance boundaries that all
three sub-commands now share. The boundary cases are deliberately exhaustive:
the off-by-one disagreements between the old per-command implementations are
precisely what this module exists to prevent.
"""

import pytest

from teloclip.overhang import (
    REASON_MAX_BREAK,
    REASON_MIN_ANCHOR,
    REASON_MIN_CLIP,
    REASON_NO_CLIP,
    AlignmentEnds,
    alignment_end,
    anchor_length,
    cigar_ops_from_string,
    cigar_ops_from_tuples,
    classify,
    classify_end,
    clip_lengths,
    ends_from_aligned_segment,
    ends_from_sam_fields,
    overhang_info_from_call,
    reference_span,
)

CONTIG_LEN = 1000


def make_ends(
    aln_start: int,
    ref_span: int,
    left_clip: int = 0,
    right_clip: int = 0,
    anchor: int = 1000,
    contig_length: int = CONTIG_LEN,
    sequence: str = '',
) -> AlignmentEnds:
    """
    Build an AlignmentEnds directly, bypassing the CIGAR adapters.

    Parameters
    ----------
    aln_start : int
        First aligned reference base, 1-based.
    ref_span : int
        Number of reference bases covered by the alignment.
    left_clip : int, optional
        Leading soft-clip length (default: 0).
    right_clip : int, optional
        Trailing soft-clip length (default: 0).
    anchor : int, optional
        Anchoring base count (default: 1000, i.e. comfortably passing).
    contig_length : int, optional
        Reference contig length (default: 1000).
    sequence : str, optional
        Read sequence (default: '').

    Returns
    -------
    AlignmentEnds
        Canonical geometry for the described alignment.
    """
    return AlignmentEnds(
        contig_name='contig1',
        contig_length=contig_length,
        aln_start=aln_start,
        aln_end=aln_start + ref_span - 1,
        left_clip=left_clip,
        right_clip=right_clip,
        anchor=anchor,
        read_name='read1',
        sequence=sequence,
    )


class TestCigarParsing:
    """CIGAR parsing helpers and the two input conventions."""

    def test_parse_simple_cigar(self):
        """A soft-clipped CIGAR splits into ordered (length, op) pairs."""
        assert cigar_ops_from_string('174M76S') == [(174, 'M'), (76, 'S')]
        assert cigar_ops_from_string('96S154M') == [(96, 'S'), (154, 'M')]

    def test_parse_multi_digit_and_all_ops(self):
        """Every SAM CIGAR operator parses, including multi-digit lengths."""
        assert cigar_ops_from_string('10S100M5I20D3N7=2X8S') == [
            (10, 'S'),
            (100, 'M'),
            (5, 'I'),
            (20, 'D'),
            (3, 'N'),
            (7, '='),
            (2, 'X'),
            (8, 'S'),
        ]

    @pytest.mark.parametrize('cigar', ['', '*'])
    def test_parse_unavailable_cigar(self, cigar):
        """An absent CIGAR yields no operations rather than raising."""
        assert cigar_ops_from_string(cigar) == []

    def test_pysam_tuples_roundtrip(self):
        """Pysam (op_code, length) tuples map onto the same representation."""
        # 4 == S, 0 == M, per 'MIDNSHP=XB'.
        assert cigar_ops_from_tuples([(4, 96), (0, 154)]) == [(96, 'S'), (154, 'M')]

    def test_pysam_tuples_none(self):
        """An unmapped read with no CIGAR yields no operations."""
        assert cigar_ops_from_tuples(None) == []

    def test_both_adapters_agree_on_span(self):
        """The raw-SAM and pysam paths compute identical geometry.

        This is the invariant that prevents the two code paths drifting apart
        the way the pre-refactor implementations did.
        """
        cigar = '10S50M5D20M8S'
        string_ops = cigar_ops_from_string(cigar)
        tuple_ops = cigar_ops_from_tuples(
            [(4, 10), (0, 50), (2, 5), (0, 20), (4, 8)]  # S M D M S
        )
        assert string_ops == tuple_ops
        assert reference_span(string_ops) == reference_span(tuple_ops)
        assert clip_lengths(string_ops) == clip_lengths(tuple_ops)


class TestSpanAndAnchor:
    """Reference span, anchor length and clip extraction."""

    def test_reference_span_counts_ref_consuming_ops_only(self):
        """M, D, N, = and X advance the reference; I, S and H do not."""
        ops = cigar_ops_from_string('10S50M5I5D3N7=2X8S')
        assert reference_span(ops) == 50 + 5 + 3 + 7 + 2

    def test_anchor_counts_matches_only(self):
        """Anchoring counts M/=/X, excluding indels and clips."""
        ops = cigar_ops_from_string('10S50M5I5D3N7=2X8S')
        assert anchor_length(ops) == 50 + 7 + 2

    def test_anchor_excludes_insertions(self):
        """Insertions are read bases with no reference support, so do not anchor.

        This is the behaviour that ``extend`` previously got wrong by computing
        the anchor as read length minus clips.
        """
        ops = cigar_ops_from_string('10S100M50I100M10S')
        assert anchor_length(ops) == 200
        # The naive read-length-minus-clips formula would have said 250.
        read_length = sum(
            length for length, op in ops if op in {'M', 'I', 'S', '=', 'X'}
        )
        assert read_length - 10 - 10 == 250

    def test_clip_lengths_reports_zero_not_none(self):
        """Absent clips report as 0 so callers can do arithmetic unguarded."""
        assert clip_lengths(cigar_ops_from_string('100M')) == (0, 0)
        assert clip_lengths(cigar_ops_from_string('50S100M')) == (50, 0)
        assert clip_lengths(cigar_ops_from_string('100M20S')) == (0, 20)
        assert clip_lengths(cigar_ops_from_string('10S100M20S')) == (10, 20)

    def test_clip_only_cigar_not_double_counted(self):
        """A single soft-clip element counts as a left clip only."""
        assert clip_lengths(cigar_ops_from_string('100S')) == (100, 0)

    def test_alignment_end_is_inclusive(self):
        """aln_end is the last aligned base, not one position past it."""
        ops = cigar_ops_from_string('100M')
        # An alignment starting at 1 spanning 100 bases ends at 100, not 101.
        assert alignment_end(1, ops) == 100

    def test_alignment_end_with_deletions(self):
        """Deletions and skips extend the reference footprint."""
        ops = cigar_ops_from_string('50M10D50M')
        assert alignment_end(1, ops) == 110


class TestGaps:
    """Gap definitions at each terminus."""

    def test_flush_alignment_has_zero_gaps(self):
        """An alignment covering the whole contig leaves no uncovered bases."""
        ends = make_ends(aln_start=1, ref_span=CONTIG_LEN)
        assert ends.gap_left == 0
        assert ends.gap_right == 0

    def test_gap_left_is_bases_before_alignment(self):
        """Starting at position 11 leaves 10 uncovered bases at the 5' end."""
        assert make_ends(aln_start=11, ref_span=100).gap_left == 10

    def test_gap_right_is_bases_after_alignment(self):
        """Ending at 990 on a 1000 bp contig leaves 10 uncovered bases."""
        assert make_ends(aln_start=891, ref_span=100).gap_right == 10


class TestMaxBreakBoundary:
    """max_break is an inclusive distance, identical on both ends."""

    @pytest.mark.parametrize('is_left', [True, False])
    def test_gap_equal_to_max_break_accepted(self, is_left):
        """Gap == max_break is inside the tolerance."""
        if is_left:
            ends = make_ends(aln_start=51, ref_span=100, left_clip=100)
        else:
            ends = make_ends(aln_start=801, ref_span=150, right_clip=100)
        call = classify_end(ends, is_left, max_break=50, min_clip=1)
        assert call.gap == 50
        assert call.accepted

    @pytest.mark.parametrize('is_left', [True, False])
    def test_gap_one_past_max_break_rejected(self, is_left):
        """Gap == max_break + 1 is outside the tolerance."""
        if is_left:
            ends = make_ends(aln_start=52, ref_span=100, left_clip=200)
        else:
            ends = make_ends(aln_start=800, ref_span=150, right_clip=200)
        call = classify_end(ends, is_left, max_break=50, min_clip=1)
        assert call.gap == 51
        assert not call.accepted
        assert call.reason == REASON_MAX_BREAK

    @pytest.mark.parametrize('is_left', [True, False])
    def test_zero_max_break_requires_flush_alignment(self, is_left):
        """max_break=0 admits only alignments reaching the terminal base."""
        if is_left:
            flush = make_ends(aln_start=1, ref_span=100, left_clip=10)
            offset = make_ends(aln_start=2, ref_span=100, left_clip=10)
        else:
            flush = make_ends(aln_start=901, ref_span=100, right_clip=10)
            offset = make_ends(aln_start=900, ref_span=100, right_clip=10)

        assert classify_end(flush, is_left, max_break=0, min_clip=1).accepted
        rejected = classify_end(offset, is_left, max_break=0, min_clip=1)
        assert not rejected.accepted
        assert rejected.reason == REASON_MAX_BREAK


class TestMinClipBoundary:
    """min_clip thresholds the overhang past the terminus, not the raw clip."""

    @pytest.mark.parametrize('is_left', [True, False])
    def test_overhang_equal_to_min_clip_accepted(self, is_left):
        """Overhang == min_clip is sufficient."""
        # gap of 10, clip of 15 => overhang of 5.
        if is_left:
            ends = make_ends(aln_start=11, ref_span=100, left_clip=15)
        else:
            ends = make_ends(aln_start=891, ref_span=100, right_clip=15)
        call = classify_end(ends, is_left, max_break=50, min_clip=5)
        assert call.overhang_len == 5
        assert call.accepted

    @pytest.mark.parametrize('is_left', [True, False])
    def test_overhang_below_min_clip_rejected(self, is_left):
        """Overhang == min_clip - 1 is insufficient."""
        if is_left:
            ends = make_ends(aln_start=11, ref_span=100, left_clip=14)
        else:
            ends = make_ends(aln_start=891, ref_span=100, right_clip=14)
        call = classify_end(ends, is_left, max_break=50, min_clip=5)
        assert call.overhang_len == 4
        assert not call.accepted
        assert call.reason == REASON_MIN_CLIP

    @pytest.mark.parametrize('is_left', [True, False])
    def test_zero_overhang_rejected(self, is_left):
        """A clip that exactly reaches the terminus adds nothing and is rejected.

        The pre-refactor right-end test accepted this case, which is how
        zero-gain extensions entered the pipeline.
        """
        if is_left:
            ends = make_ends(aln_start=11, ref_span=100, left_clip=10)
        else:
            ends = make_ends(aln_start=891, ref_span=100, right_clip=10)
        call = classify_end(ends, is_left, max_break=50, min_clip=1)
        assert call.overhang_len == 0
        assert not call.accepted
        assert call.reason == REASON_MIN_CLIP

    @pytest.mark.parametrize('is_left', [True, False])
    def test_negative_overhang_rejected(self, is_left):
        """A clip stopping short of the terminus would shorten the contig.

        gap 20, clip 3: applying this would trim 20 bases and add 3, for a net
        loss of 17. This is the exact shape of the bug this module removes.
        """
        if is_left:
            ends = make_ends(aln_start=21, ref_span=100, left_clip=3)
        else:
            ends = make_ends(aln_start=881, ref_span=100, right_clip=3)
        call = classify_end(ends, is_left, max_break=50, min_clip=1)
        assert call.overhang_len == -17
        assert not call.accepted
        assert call.reason == REASON_MIN_CLIP


class TestNoClip:
    """Ends with no soft clip are rejected before any threshold applies."""

    @pytest.mark.parametrize('is_left', [True, False])
    def test_absent_clip_rejected(self, is_left):
        """An unclipped end is reported as no_clip, not as min_clip."""
        ends = make_ends(aln_start=1, ref_span=CONTIG_LEN)
        call = classify_end(ends, is_left, max_break=50, min_clip=1)
        assert not call.accepted
        assert call.reason == REASON_NO_CLIP


class TestSymmetry:
    """The two ends must apply identical rules."""

    @pytest.mark.parametrize(
        'gap,clip', [(0, 10), (5, 10), (10, 10), (10, 11), (50, 60), (51, 60), (20, 3)]
    )
    def test_mirrored_cases_give_identical_verdicts(self, gap, clip):
        """A case and its mirror image are judged the same way.

        Any left/right asymmetry in the acceptance rules shows up here.
        """
        left = make_ends(aln_start=gap + 1, ref_span=100, left_clip=clip)
        right = make_ends(
            aln_start=CONTIG_LEN - gap - 100 + 1, ref_span=100, right_clip=clip
        )

        left_call = classify_end(left, True, max_break=50, min_clip=1)
        right_call = classify_end(right, False, max_break=50, min_clip=1)

        assert left_call.gap == right_call.gap == gap
        assert left_call.overhang_len == right_call.overhang_len
        assert left_call.accepted == right_call.accepted
        assert left_call.reason == right_call.reason


class TestNetGainInvariant:
    """Every accepted call must represent a real increase in contig length."""

    @pytest.mark.parametrize('is_left', [True, False])
    def test_accepted_implies_positive_net_gain(self, is_left):
        """Accepted => net_gain >= min_clip >= 1, so the contig never shrinks.

        Swept across the full grid of gaps, clip lengths and min_clip settings
        rather than parametrised, to keep one property to one test.
        """
        checked_accepted = 0

        for min_clip in (1, 2, 5, 50):
            for gap in range(0, 12):
                for clip in range(0, 60, 7):
                    if is_left:
                        ends = make_ends(
                            aln_start=gap + 1, ref_span=100, left_clip=clip
                        )
                    else:
                        ends = make_ends(
                            aln_start=CONTIG_LEN - gap - 100 + 1,
                            ref_span=100,
                            right_clip=clip,
                        )
                    call = classify_end(ends, is_left, max_break=10, min_clip=min_clip)
                    if call.accepted:
                        checked_accepted += 1
                        assert call.net_gain >= min_clip >= 1, (
                            f'min_clip={min_clip} gap={gap} clip={clip}'
                        )

        # Guard against the sweep silently accepting nothing and passing vacuously.
        assert checked_accepted > 0

    def test_trim_len_never_negative(self):
        """An alignment overhanging the contig end reports no trim."""
        # ref_span runs past the contig end, so gap_right is negative.
        ends = make_ends(aln_start=951, ref_span=100, right_clip=20)
        call = classify_end(ends, False, max_break=50, min_clip=1)
        assert call.gap == -50
        assert call.trim_len == 0


class TestOverhangSequence:
    """Clip sequence extraction."""

    def test_left_clip_taken_from_read_start(self):
        """The left overhang is the leading clip, in full."""
        ends = make_ends(aln_start=1, ref_span=6, left_clip=4, sequence='AAAATTTTTT')
        call = classify_end(ends, True, max_break=50, min_clip=1)
        assert call.overhang_sequence(ends) == 'AAAA'

    def test_right_clip_taken_from_read_end(self):
        """The right overhang is the trailing clip, in full."""
        ends = make_ends(
            aln_start=CONTIG_LEN - 5, ref_span=6, right_clip=4, sequence='TTTTTTGGGG'
        )
        call = classify_end(ends, False, max_break=50, min_clip=1)
        assert call.overhang_sequence(ends) == 'GGGG'

    def test_full_clip_returned_not_just_novel_tail(self):
        """The whole clip is grafted; the trimmed gap is replaced by read bases.

        Net length change still equals net_gain: 4 bases added, 2 trimmed.
        """
        ends = make_ends(aln_start=3, ref_span=6, left_clip=4, sequence='AAAATTTTTT')
        call = classify_end(ends, True, max_break=50, min_clip=1)
        assert call.overhang_sequence(ends) == 'AAAA'
        assert call.trim_len == 2
        assert call.net_gain == 2

    def test_missing_sequence_yields_empty(self):
        """A record without SEQ yields an empty overhang rather than raising."""
        ends = make_ends(aln_start=1, ref_span=100, left_clip=10)
        call = classify_end(ends, True, max_break=50, min_clip=1)
        assert call.overhang_sequence(ends) == ''


class TestAdapters:
    """Building AlignmentEnds from each input representation."""

    def test_from_sam_fields(self):
        """A tab-split SAM record maps onto canonical geometry."""
        fields = [
            'read1',
            '0',
            'contig1',
            '11',
            '60',
            '10S100M20S',
            '*',
            '0',
            '0',
            'A' * 130,
            'I' * 130,
        ]
        ends = ends_from_sam_fields(fields, CONTIG_LEN)

        assert ends.read_name == 'read1'
        assert ends.contig_name == 'contig1'
        assert ends.aln_start == 11
        assert ends.aln_end == 110  # inclusive: 11 + 100 - 1
        assert ends.left_clip == 10
        assert ends.right_clip == 20
        assert ends.anchor == 100
        assert ends.gap_left == 10
        assert ends.gap_right == CONTIG_LEN - 110

    def test_from_aligned_segment(self):
        """A pysam-like segment maps onto identical geometry."""

        class FakeSegment:
            """Minimal stand-in exposing the attributes the adapter reads."""

            reference_start = 10  # 0-based; equivalent to POS 11
            cigartuples = [(4, 10), (0, 100), (4, 20)]  # 10S100M20S
            query_name = 'read1'
            query_sequence = 'A' * 130

        ends = ends_from_aligned_segment(FakeSegment(), CONTIG_LEN, 'contig1')

        assert ends.aln_start == 11
        assert ends.aln_end == 110
        assert ends.left_clip == 10
        assert ends.right_clip == 20
        assert ends.anchor == 100

    def test_adapters_produce_equal_ends(self):
        """The two adapters agree on every field for the same alignment."""
        fields = [
            'read1',
            '0',
            'contig1',
            '11',
            '60',
            '10S100M20S',
            '*',
            '0',
            '0',
            'A' * 130,
            'I' * 130,
        ]

        class FakeSegment:
            """Same alignment, expressed the way pysam would report it."""

            reference_start = 10
            cigartuples = [(4, 10), (0, 100), (4, 20)]
            query_name = 'read1'
            query_sequence = 'A' * 130

        assert ends_from_sam_fields(fields, CONTIG_LEN) == ends_from_aligned_segment(
            FakeSegment(), CONTIG_LEN, 'contig1'
        )


class TestClassifyBothEnds:
    """The two-end entry point, including the anchor gate."""

    def test_both_ends_classified(self):
        """A read clipped at both ends of a short contig yields two calls."""
        ends = make_ends(aln_start=1, ref_span=CONTIG_LEN, left_clip=50, right_clip=60)
        left, right = classify(ends, max_break=10, min_clip=1, min_anchor=0)
        assert left.accepted and left.net_gain == 50
        assert right.accepted and right.net_gain == 60

    def test_min_anchor_rejects_both_ends(self):
        """A poorly anchored read is rejected at both ends with one reason."""
        ends = make_ends(
            aln_start=1, ref_span=CONTIG_LEN, left_clip=50, right_clip=60, anchor=99
        )
        left, right = classify(ends, max_break=10, min_clip=1, min_anchor=100)
        assert not left.accepted and left.reason == REASON_MIN_ANCHOR
        assert not right.accepted and right.reason == REASON_MIN_ANCHOR

    def test_min_anchor_boundary_is_inclusive(self):
        """Anchor == min_anchor passes the gate."""
        ends = make_ends(aln_start=1, ref_span=CONTIG_LEN, left_clip=50, anchor=100)
        left, _ = classify(ends, max_break=10, min_clip=1, min_anchor=100)
        assert left.accepted


class TestOverhangInfoAdapter:
    """Conversion to the OverhangInfo record consumed by the extend path."""

    def test_populates_canonical_fields(self):
        """The record carries inclusive coordinates and the pre-computed gain."""
        ends = make_ends(
            aln_start=3, ref_span=6, left_clip=4, anchor=6, sequence='AAAATTTTTT'
        )
        call = classify_end(ends, True, max_break=50, min_clip=1)
        info = overhang_info_from_call(ends, call)

        assert info.sequence == 'AAAA'
        assert info.length == 4
        assert info.clip_length == 4
        assert info.alignment_pos == 3
        assert info.alignment_end == 8  # inclusive
        assert info.anchor_length == 6
        assert info.contig_name == 'contig1'
        assert info.read_name == 'read1'
        assert info.is_left is True
        assert info.net_gain == 2
