"""
Canonical definition of a soft-clip overhang at a contig terminus.

This module is the single source of truth shared by ``teloclip filter``,
``teloclip extract`` and ``teloclip extend``. Historically each sub-command
carried its own copy of the "is this clip a valid terminal overhang?" test, and
the three copies disagreed with one another by one base on both ends and about
whether ``min_clip`` applied at all. Everything now routes through
:func:`classify`.

Notes
-----
**Coordinate convention.**
All reference coordinates here are **1-based and inclusive**. This convention is
required by :func:`teloclip.core.extension.calculate_extension_position`, which
derives the right-hand trim as ``contig_length - alignment_end``; that is only
correct when ``alignment_end`` is the last aligned reference base rather than one
position past it.

For an alignment starting at reference position ``p`` on a contig of length
``L``::

    ref_span  = sum of CIGAR lengths over {M, D, N, =, X}
    aln_start = p                     # first aligned reference base
    aln_end   = p + ref_span - 1      # last aligned reference base, inclusive

    gap_left  = aln_start - 1         # contig bases at the 5' terminus not covered
    gap_right = L - aln_end           # contig bases at the 3' terminus not covered

    overhang  = clip_len - gap        # read bases lying *past* the terminus

The ``overhang`` formula follows from the soft clip being contiguous with the
alignment in read space: the clipped bases notionally occupy reference positions
``aln_start - 1``, ``aln_start - 2``, ... ``aln_start - clip_len``. Of those,
``gap_left`` fall inside the contig and the remainder fall outside it. The right
end is the mirror image.

**Acceptance criteria.**
Applied identically to both ends::

    max_break:  gap      <= max_break   # inclusive; 0 means flush with the terminus
    min_clip:   overhang >= min_clip    # measured past the terminus, not raw clip length
    min_anchor: anchor   >= min_anchor  # anchor counts M/=/X bases only

``min_clip`` deliberately thresholds the *overhang* rather than the raw clip
length, matching the long-standing help text ("require clip to extend past ref
contig end by at least N bases"). A 500 bp clip sitting 499 bp inside the contig
contributes exactly 1 bp of novel sequence, and that is the biologically
meaningful quantity.

Because callers clamp ``min_clip`` to at least 1, every accepted call carries
``net_gain >= 1``. An extension that would shorten a contig is therefore not
representable, rather than merely unlikely.
"""

from dataclasses import dataclass
import re
from typing import Iterable, List, Optional, Sequence, Tuple

# CIGAR operations that advance the position on the reference sequence.
REF_CONSUMING = frozenset({'M', 'D', 'N', '=', 'X'})

# CIGAR operations that count as anchoring evidence. Deletions and skips advance
# the reference but contribute no read bases, and insertions contribute read
# bases with no reference support, so none of them anchor the read.
ANCHOR_OPS = frozenset({'M', '=', 'X'})

# pysam encodes CIGAR operations as integers indexing this string.
_PYSAM_CIGAR_CODES = 'MIDNSHP=XB'

# Rejection reasons. These strings double as exclusion-statistic keys so that
# every sub-command reports the same buckets.
REASON_NO_CLIP = 'no_clip'
REASON_MAX_BREAK = 'max_break'
REASON_MIN_CLIP = 'min_clip'
REASON_MIN_ANCHOR = 'min_anchor'

# A CIGAR parsed into (length, operation) pairs, e.g. [(174, 'M'), (76, 'S')].
CigarOps = List[Tuple[int, str]]

_CIGAR_ELEMENT = re.compile(r'([0-9]+)([MIDNSHP=X])')


def cigar_ops_from_string(cigar: str) -> CigarOps:
    """
    Parse a CIGAR string into (length, operation) pairs.

    Parameters
    ----------
    cigar : str
        CIGAR string from a SAM record, e.g. ``'96S154M'``.

    Returns
    -------
    CigarOps
        List of ``(length, operation)`` tuples in CIGAR order. An empty or
        unavailable CIGAR (``''`` or ``'*'``) yields an empty list.

    Examples
    --------
    >>> cigar_ops_from_string('174M76S')
    [(174, 'M'), (76, 'S')]
    """
    if not cigar or cigar == '*':
        return []
    return [(int(length), op) for length, op in _CIGAR_ELEMENT.findall(cigar)]


def cigar_ops_from_tuples(cigartuples: Optional[Iterable[Tuple[int, int]]]) -> CigarOps:
    """
    Convert pysam-style CIGAR tuples into (length, operation) pairs.

    pysam reports CIGAR elements as ``(operation_code, length)`` where the code
    indexes ``'MIDNSHP=XB'``. This is the inverse mapping, so that the pysam and
    raw-SAM code paths share every downstream calculation.

    Parameters
    ----------
    cigartuples : iterable of tuple of int, optional
        Sequence of ``(operation_code, length)`` pairs, as returned by
        ``pysam.AlignedSegment.cigartuples``. ``None`` yields an empty list.

    Returns
    -------
    CigarOps
        List of ``(length, operation)`` tuples in CIGAR order.

    Examples
    --------
    >>> cigar_ops_from_tuples([(4, 96), (0, 154)])
    [(96, 'S'), (154, 'M')]
    """
    if not cigartuples:
        return []
    return [
        (length, _PYSAM_CIGAR_CODES[op])
        for op, length in cigartuples
        if 0 <= op < len(_PYSAM_CIGAR_CODES)
    ]


def reference_span(ops: CigarOps) -> int:
    """
    Calculate how many reference bases an alignment covers.

    Parameters
    ----------
    ops : CigarOps
        Parsed CIGAR operations.

    Returns
    -------
    int
        Sum of the lengths of all reference-consuming operations
        (``M``, ``D``, ``N``, ``=``, ``X``).
    """
    return sum(length for length, op in ops if op in REF_CONSUMING)


def anchor_length(ops: CigarOps) -> int:
    """
    Calculate the number of read bases anchored to the reference.

    Only ``M``, ``=`` and ``X`` count. Insertions are read bases with no
    reference support and deletions are reference bases with no read support, so
    neither provides evidence that the read is correctly placed.

    Parameters
    ----------
    ops : CigarOps
        Parsed CIGAR operations.

    Returns
    -------
    int
        Number of anchoring bases.
    """
    return sum(length for length, op in ops if op in ANCHOR_OPS)


def clip_lengths(ops: CigarOps) -> Tuple[int, int]:
    """
    Get the leading and trailing soft-clip lengths of an alignment.

    Parameters
    ----------
    ops : CigarOps
        Parsed CIGAR operations.

    Returns
    -------
    tuple of int
        ``(left_clip, right_clip)``. Unlike
        :func:`teloclip.io.sam.checkClips`, absent clips are reported as ``0``
        rather than ``None``, so callers can do arithmetic without guarding.

    Examples
    --------
    >>> clip_lengths([(10, 'S'), (100, 'M'), (20, 'S')])
    (10, 20)
    >>> clip_lengths([(100, 'M')])
    (0, 0)
    """
    if not ops:
        return (0, 0)
    left = ops[0][0] if ops[0][1] == 'S' else 0
    # A single-element soft-clip-only CIGAR must not be counted at both ends.
    right = ops[-1][0] if len(ops) > 1 and ops[-1][1] == 'S' else 0
    return (left, right)


def alignment_end(aln_start: int, ops: CigarOps) -> int:
    """
    Calculate the last aligned reference position, 1-based and inclusive.

    Parameters
    ----------
    aln_start : int
        First aligned reference position (1-based).
    ops : CigarOps
        Parsed CIGAR operations.

    Returns
    -------
    int
        ``aln_start + reference_span(ops) - 1``.
    """
    return aln_start + reference_span(ops) - 1


@dataclass(frozen=True)
class AlignmentEnds:
    """
    Terminus-relevant geometry of a single alignment, in canonical coordinates.

    Built by one of the adapter functions (:func:`ends_from_sam_fields` or
    :func:`ends_from_aligned_segment`) so that the raw-SAM and pysam code paths
    converge before any acceptance decision is made.
    """

    contig_name: str
    contig_length: int
    aln_start: int
    aln_end: int
    left_clip: int
    right_clip: int
    anchor: int
    read_name: str
    sequence: str = ''
    # Retained for display and for CIGAR-aware layout in the HTML report.
    # None of these take part in the acceptance decision.
    cigar: str = ''
    mapq: int = -1
    flag: int = 0

    @property
    def gap_left(self) -> int:
        """
        Contig bases at the 5' terminus not covered by this alignment.

        Returns
        -------
        int
            ``aln_start - 1``.
        """
        return self.aln_start - 1

    @property
    def gap_right(self) -> int:
        """
        Contig bases at the 3' terminus not covered by this alignment.

        Returns
        -------
        int
            ``contig_length - aln_end``.
        """
        return self.contig_length - self.aln_end


@dataclass(frozen=True)
class OverhangCall:
    """
    The verdict for one end of one alignment.

    A call is produced for both ends of every alignment, whether accepted or
    not, so that rejection reasons can be tallied consistently.
    """

    is_left: bool
    gap: int
    clip_len: int
    overhang_len: int
    accepted: bool
    reason: Optional[str] = None

    @property
    def net_gain(self) -> int:
        """
        Change in contig length if this overhang were applied.

        Applying an overhang trims ``gap`` bases from the terminus and grafts on
        the full ``clip_len``-base clip, so the net change is exactly the number
        of clipped bases lying past the terminus.

        Returns
        -------
        int
            ``overhang_len``. Negative or zero for a clip that does not reach
            past the contig terminus.
        """
        return self.overhang_len

    @property
    def trim_len(self) -> int:
        """
        Contig bases removed at this terminus if the overhang were applied.

        Returns
        -------
        int
            ``max(0, gap)``.
        """
        return max(0, self.gap)

    def overhang_sequence(self, ends: AlignmentEnds) -> str:
        """
        Extract the clip sequence to graft onto the contig.

        The *entire* soft clip is returned, not just the ``overhang_len`` bases
        past the terminus: the ``gap`` contig bases are trimmed away and
        replaced by the read's own version of that region, plus the novel tail.
        The net length change remains :attr:`net_gain`.

        Parameters
        ----------
        ends : AlignmentEnds
            The alignment this call was derived from.

        Returns
        -------
        str
            Clip sequence of length ``clip_len``, or ``''`` if the read
            sequence was unavailable.
        """
        if not ends.sequence or self.clip_len <= 0:
            return ''
        if self.is_left:
            return ends.sequence[: self.clip_len]
        return ends.sequence[-self.clip_len :]

    def read_slice(self, ends: AlignmentEnds, limit: int) -> Tuple[str, int]:
        """
        Extract the portion of the read needed to render this terminus.

        The HTML report walks the CIGAR to align the read against the contig,
        which needs real read bases rather than a pre-extracted anchor. Keeping
        the whole read would be wasteful for a 100 kb long read, so only the
        end abutting the terminus is kept, along with how many bases were
        dropped from the front so CIGAR offsets still resolve.

        Parameters
        ----------
        ends : AlignmentEnds
            The alignment this call was derived from.
        limit : int
            Maximum read bases to retain. Zero or less returns nothing, which
            is the default for runs that do not need it.

        Returns
        -------
        tuple
            ``(sequence, offset)`` where ``offset`` is the number of read bases
            dropped from the start of the record.
        """
        if limit <= 0 or not ends.sequence:
            return '', 0
        if len(ends.sequence) <= limit:
            return ends.sequence, 0
        if self.is_left:
            # The 5' terminus is served by the start of the read.
            return ends.sequence[:limit], 0
        # The 3' terminus is served by the end of it.
        offset = len(ends.sequence) - limit
        return ends.sequence[offset:], offset


# SAM field indices, per the SAM specification.
_SAM_QNAME = 0
_SAM_FLAG = 1
_SAM_RNAME = 2
_SAM_POS = 3
_SAM_MAPQ = 4
_SAM_CIGAR = 5
_SAM_SEQ = 9


def ends_from_sam_fields(fields: Sequence[str], contig_length: int) -> AlignmentEnds:
    """
    Build :class:`AlignmentEnds` from a tab-split SAM line.

    Used by the ``filter`` and ``extract`` code paths, which stream raw SAM text
    rather than going through pysam.

    Parameters
    ----------
    fields : sequence of str
        A SAM record already split on tabs. At least 10 fields are required.
    contig_length : int
        Length of the reference contig this record aligns to.

    Returns
    -------
    AlignmentEnds
        Canonical geometry for this alignment.

    Raises
    ------
    IndexError
        If ``fields`` is too short to contain a SEQ column.
    ValueError
        If the POS column is not an integer.
    """
    ops = cigar_ops_from_string(fields[_SAM_CIGAR])
    aln_start = int(fields[_SAM_POS])
    left_clip, right_clip = clip_lengths(ops)

    return AlignmentEnds(
        contig_name=fields[_SAM_RNAME],
        contig_length=contig_length,
        aln_start=aln_start,
        aln_end=alignment_end(aln_start, ops),
        left_clip=left_clip,
        right_clip=right_clip,
        anchor=anchor_length(ops),
        read_name=fields[_SAM_QNAME],
        sequence=fields[_SAM_SEQ],
        cigar=fields[_SAM_CIGAR],
        mapq=int(fields[_SAM_MAPQ]) if fields[_SAM_MAPQ].isdigit() else -1,
        flag=int(fields[_SAM_FLAG]),
    )


def ends_from_aligned_segment(
    aln, contig_length: int, contig_name: str
) -> AlignmentEnds:
    """
    Build :class:`AlignmentEnds` from a pysam ``AlignedSegment``.

    Used by the ``extend`` code path. The segment is duck-typed on
    ``reference_start``, ``cigartuples``, ``query_name`` and ``query_sequence``,
    so this module never imports pysam.

    ``reference_end`` is deliberately not used: computing the span through
    :func:`reference_span` on both code paths is what guarantees the raw-SAM and
    pysam routes can never drift apart.

    Parameters
    ----------
    aln : object
        A ``pysam.AlignedSegment`` (or anything exposing the same attributes).
    contig_length : int
        Length of the reference contig this record aligns to.
    contig_name : str
        Name of the reference contig.

    Returns
    -------
    AlignmentEnds
        Canonical geometry for this alignment.
    """
    ops = cigar_ops_from_tuples(aln.cigartuples)
    # pysam reports a 0-based start; canonical coordinates are 1-based.
    aln_start = aln.reference_start + 1
    left_clip, right_clip = clip_lengths(ops)

    return AlignmentEnds(
        contig_name=contig_name,
        contig_length=contig_length,
        aln_start=aln_start,
        aln_end=alignment_end(aln_start, ops),
        left_clip=left_clip,
        right_clip=right_clip,
        anchor=anchor_length(ops),
        read_name=aln.query_name or '',
        sequence=aln.query_sequence or '',
        cigar=aln.cigarstring or '',
        mapq=-1 if aln.mapping_quality is None else int(aln.mapping_quality),
        flag=int(aln.flag),
    )


def classify_end(
    ends: AlignmentEnds,
    is_left: bool,
    *,
    max_break: int,
    min_clip: int,
) -> OverhangCall:
    """
    Decide whether one end of one alignment is a valid terminal overhang.

    The end is accepted when all of the following hold::

        clip_len > 0
        gap <= max_break
        (clip_len - gap) >= min_clip

    Checks are evaluated in the fixed order no-clip, then ``max_break``, then
    ``min_clip``, so that rejection-reason tallies are deterministic and
    directly comparable between sub-commands.

    Parameters
    ----------
    ends : AlignmentEnds
        Canonical geometry for the alignment.
    is_left : bool
        ``True`` to test the 5' terminus, ``False`` for the 3' terminus.
    max_break : int
        Maximum tolerated number of uncovered contig bases between the terminus
        and the start of the alignment. Inclusive: ``gap == max_break`` passes.
    min_clip : int
        Minimum number of clipped bases required to lie past the terminus.

    Returns
    -------
    OverhangCall
        The verdict, carrying the geometry needed to apply or tally it.
    """
    if is_left:
        gap = ends.gap_left
        clip_len = ends.left_clip
    else:
        gap = ends.gap_right
        clip_len = ends.right_clip

    overhang_len = clip_len - gap

    if clip_len <= 0:
        reason = REASON_NO_CLIP
    elif gap > max_break:
        reason = REASON_MAX_BREAK
    elif overhang_len < min_clip:
        reason = REASON_MIN_CLIP
    else:
        reason = None

    return OverhangCall(
        is_left=is_left,
        gap=gap,
        clip_len=clip_len,
        overhang_len=overhang_len,
        accepted=reason is None,
        reason=reason,
    )


def classify(
    ends: AlignmentEnds,
    *,
    max_break: int,
    min_clip: int,
    min_anchor: int = 0,
) -> Tuple[OverhangCall, OverhangCall]:
    """
    Classify both ends of an alignment against the terminal-overhang criteria.

    Parameters
    ----------
    ends : AlignmentEnds
        Canonical geometry for the alignment.
    max_break : int
        Maximum tolerated gap between a terminus and the alignment (inclusive).
    min_clip : int
        Minimum number of clipped bases required past the terminus.
    min_anchor : int, optional
        Minimum number of anchoring (``M``/``=``/``X``) bases required. If the
        alignment falls short, both ends are rejected with
        :data:`REASON_MIN_ANCHOR` and no other check is applied (default: 0).

    Returns
    -------
    tuple of OverhangCall
        ``(left_call, right_call)``.
    """
    if ends.anchor < min_anchor:
        return (
            OverhangCall(
                is_left=True,
                gap=ends.gap_left,
                clip_len=ends.left_clip,
                overhang_len=ends.left_clip - ends.gap_left,
                accepted=False,
                reason=REASON_MIN_ANCHOR,
            ),
            OverhangCall(
                is_left=False,
                gap=ends.gap_right,
                clip_len=ends.right_clip,
                overhang_len=ends.right_clip - ends.gap_right,
                accepted=False,
                reason=REASON_MIN_ANCHOR,
            ),
        )

    return (
        classify_end(ends, True, max_break=max_break, min_clip=min_clip),
        classify_end(ends, False, max_break=max_break, min_clip=min_clip),
    )


def overhang_info_from_call(
    ends: AlignmentEnds, call: OverhangCall, anchor_context: int = 0
) -> 'OverhangInfo':  # noqa: F821
    """
    Adapt an accepted call to the :class:`teloclip.core.analysis.OverhangInfo` record.

    Parameters
    ----------
    ends : AlignmentEnds
        Canonical geometry for the alignment.
    call : OverhangCall
        An accepted call for one end of that alignment.
    anchor_context : int, optional
        Aligned bases adjacent to the clip to retain for later display
        (default: 0, i.e. none). Only the HTML report needs these.

    Returns
    -------
    OverhangInfo
        Populated overhang record, carrying the canonical inclusive
        ``alignment_end`` and the pre-computed ``net_gain``.
    """
    # Imported lazily: analysis imports samops, which imports this module.
    from .analysis import OverhangInfo

    _slice = call.read_slice(ends, anchor_context)

    return OverhangInfo(
        sequence=call.overhang_sequence(ends),
        length=call.clip_len,
        alignment_pos=ends.aln_start,
        alignment_end=ends.aln_end,
        read_name=ends.read_name,
        is_left=call.is_left,
        clip_length=call.clip_len,
        anchor_length=ends.anchor,
        contig_name=ends.contig_name,
        net_gain=call.net_gain,
        read_seq=_slice[0],
        read_seq_offset=_slice[1],
        cigar=ends.cigar,
        mapq=ends.mapq,
        flag=ends.flag,
    )
