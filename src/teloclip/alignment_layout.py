"""
Column layout for the overhang alignment view.

Laying reads out against a contig terminus is a small multiple-sequence
alignment problem, and doing it naively gets it wrong. Two things have to be
handled:

**The CIGAR.** A read's anchored portion is not a flat copy of the reference.
Deletions mean the read has no base where the reference does; insertions mean
the read has bases the reference does not. Pasting the read sequence next to the
reference and hoping they line up drifts by one column per indel, which makes a
perfectly good alignment look like a mismatched one.

**Insertions from other reads.** An insertion is a column that exists in one
read and in no other. Once any read introduces one, every other row — including
the reference — needs a gap at that column, or the rows below stop lining up.
So the column set has to be computed across all reads at an end before any row
can be rendered.

Coordinates are offsets from the contig terminus: 0 is the terminal base,
negative is outside the contig at the left end, positive outside it at the
right. Within one reference offset, insertion columns are ordered after the
reference base itself.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .overhang import cigar_ops_from_string

# Column identity. ``ins`` is 0 for a real reference column and 1..n for the
# insertion columns that sit after it.
Column = Tuple[int, int]

# What a row holds at one column.
GAP = '-'
BLANK = ' '


@dataclass
class Placed:
    """A read's bases keyed by the column each occupies."""

    #: Column to base, for the anchored (CIGAR-aligned) portion.
    anchor: Dict[Column, str] = field(default_factory=dict)
    #: Column to base, for the soft-clipped portion.
    clip: Dict[Column, str] = field(default_factory=dict)
    #: Reference offsets the alignment covers, for drawing a continuous row.
    first_col: Optional[Column] = None
    last_col: Optional[Column] = None


def _iter_anchor(
    cigar_ops: Sequence[Tuple[int, str]],
    sequence: str,
    aln_start: int,
    terminus: int,
    lo: int,
    hi: int,
    seq_offset: int = 0,
):
    """
    Walk a CIGAR, yielding each anchored base with its column.

    Match runs are clipped to the rendered window before iterating, so a read
    with a 100 kb alignment costs the same as one with a 500 bp alignment.

    Parameters
    ----------
    cigar_ops : sequence of tuple
        Parsed ``(length, operation)`` pairs for the whole read.
    sequence : str
        Full read sequence as it appears in the SAM record.
    aln_start : int
        First aligned reference position, 1-based.
    terminus : int
        Reference position that offset 0 corresponds to.
    lo : int
        First reference position of the window, inclusive.
    hi : int
        Last reference position of the window, inclusive.
    seq_offset : int, optional
        Read bases dropped from the front of ``sequence`` (default: 0). CIGAR
        positions are expressed against the full record, so they are shifted by
        this before indexing.

    Yields
    ------
    tuple
        ``(ref_offset, ins_index, base)``. ``ins_index`` is 0 for a base
        aligned to a reference position and 1..n for inserted bases, which sit
        after the position they follow. Deletions yield nothing, leaving those
        columns empty for this read.
    """
    read_i = 0
    ref_pos = aln_start
    end = seq_offset + len(sequence)

    def base(index: int) -> str:
        """
        Read one base by its index in the full record.

        Parameters
        ----------
        index : int
            Position in the untruncated read.

        Returns
        -------
        str
            The base, or '' when it falls outside the retained slice.
        """
        return sequence[index - seq_offset] if seq_offset <= index < end else ''

    for length, op in cigar_ops:
        if op == 'S':
            # Soft clip consumes read, not reference; handled separately.
            read_i += length
        elif op == 'H':
            # Hard clip consumes neither; those bases are not in the record.
            continue
        elif op == 'I':
            # Insertion: read bases with no reference position. They attach to
            # the reference position just consumed.
            anchor_at = ref_pos - 1
            if lo <= anchor_at <= hi:
                for k in range(length):
                    b = base(read_i + k)
                    if b:
                        yield (anchor_at - terminus, k + 1, b)
            read_i += length
        elif op in ('D', 'N'):
            # Reference advances, read does not.
            ref_pos += length
        elif op in ('M', '=', 'X'):
            # Clip the run to the window rather than walking all of it.
            start = max(ref_pos, lo)
            stop = min(ref_pos + length - 1, hi)
            for pos in range(start, stop + 1):
                b = base(read_i + (pos - ref_pos))
                if b:
                    yield (pos - terminus, 0, b)
            read_i += length
            ref_pos += length
            # No early exit here: an insertion immediately following the last
            # in-window match still attaches to a rendered column, and would be
            # lost. The per-base loop above is already clipped to the window, so
            # walking the rest of the CIGAR costs nothing.
        # P and anything unrecognised consume neither read nor reference.


def place_read(
    *,
    cigar: str,
    sequence: str,
    aln_start: int,
    aln_end: int,
    contig_length: int,
    clip_len: int,
    is_left: bool,
    window: int,
    max_overhang: int,
    seq_offset: int = 0,
) -> Placed:
    """
    Place one read's bases into columns relative to the contig terminus.

    Parameters
    ----------
    cigar : str
        CIGAR string from the SAM record.
    sequence : str
        Full read sequence from the SAM record.
    aln_start : int
        First aligned reference position, 1-based.
    aln_end : int
        Last aligned reference position, 1-based inclusive.
    contig_length : int
        Length of the contig.
    clip_len : int
        Length of the soft clip at the terminus being rendered.
    is_left : bool
        True for a 5'-terminus overhang.
    window : int
        Reference positions of anchor to render.
    max_overhang : int
        Maximum clipped bases to render.
    seq_offset : int, optional
        Read bases dropped from the front of ``sequence`` (default: 0).

    Returns
    -------
    Placed
        The read's anchored and clipped bases, keyed by column.
    """
    ops = cigar_ops_from_string(cigar)
    placed = Placed()

    terminus = 1 if is_left else contig_length

    if not sequence:
        return placed

    if is_left:
        lo, hi = terminus, terminus + window - 1
    else:
        lo, hi = terminus - window + 1, terminus

    for offset, ins, base in _iter_anchor(
        ops, sequence, aln_start, terminus, lo, hi, seq_offset
    ):
        placed.anchor[(offset, ins)] = base

    # The soft clip is unaligned, so it has no CIGAR structure of its own: its
    # bases simply run outward from where the alignment starts or ends.
    if clip_len > 0:
        if is_left:
            clip_seq = sequence[: max(0, clip_len - seq_offset)]
            if max_overhang and len(clip_seq) > max_overhang:
                clip_seq = clip_seq[-max_overhang:]
            # Last clip base sits immediately before the first aligned base.
            end_offset = (aln_start - terminus) - 1
            start_offset = end_offset - len(clip_seq) + 1
            for k, base in enumerate(clip_seq):
                placed.clip[(start_offset + k, 0)] = base
        else:
            clip_seq = sequence[-clip_len:] if clip_len <= len(sequence) else sequence
            if max_overhang and len(clip_seq) > max_overhang:
                clip_seq = clip_seq[:max_overhang]
            # First clip base sits immediately after the last aligned base.
            start_offset = (aln_end - terminus) + 1
            for k, base in enumerate(clip_seq):
                placed.clip[(start_offset + k, 0)] = base

    occupied = list(placed.anchor) + list(placed.clip)
    if occupied:
        placed.first_col = min(occupied)
        placed.last_col = max(occupied)

    return placed


def build_columns(
    placements: Sequence[Placed],
    reference_offsets: Sequence[int],
) -> List[Column]:
    """
    Compute the shared column order for a set of placed reads.

    Every reference offset in the window gets a column, plus one column for
    each insertion any read introduced at that offset. Because the set is
    computed across all reads before rendering, a gap opened by one read is
    present in every row, which is what keeps the rows aligned.

    Parameters
    ----------
    placements : sequence of Placed
        Reads already placed by :func:`place_read`.
    reference_offsets : sequence of int
        Offsets the reference row occupies.

    Returns
    -------
    list of Column
        Columns in render order.
    """
    # Widest insertion seen at each reference offset.
    max_ins: Dict[int, int] = {}
    offsets = set(reference_offsets)

    for placed in placements:
        for offset, ins in list(placed.anchor) + list(placed.clip):
            offsets.add(offset)
            if ins > max_ins.get(offset, 0):
                max_ins[offset] = ins

    columns: List[Column] = []
    for offset in sorted(offsets):
        columns.append((offset, 0))
        for k in range(1, max_ins.get(offset, 0) + 1):
            columns.append((offset, k))
    return columns


def render_row(
    placed: Placed,
    columns: Sequence[Column],
) -> Tuple[int, List[Tuple[str, str]]]:
    """
    Render one read as a run of styled segments over the shared columns.

    Clip and anchor are returned as separate runs rather than one string, so
    the caller can colour them differently without re-deriving where the
    boundary fell once gaps are in play.

    Parameters
    ----------
    placed : Placed
        The read's placed bases.
    columns : sequence of Column
        Shared column order from :func:`build_columns`.

    Returns
    -------
    tuple
        ``(lead, runs)`` where ``lead`` is the number of blank columns before
        the read begins and ``runs`` is a list of ``(kind, text)`` pairs with
        kind in ``{'clip', 'anchor', 'gap'}``. A ``gap`` run is a column the
        read does not cover: a deletion in the read, or an insertion belonging
        to another read.
    """
    if placed.first_col is None:
        return 0, []

    lead = 0
    started = False
    kinds: List[str] = []
    chars: List[str] = []

    for col in columns:
        if not started:
            if col < placed.first_col:
                lead += 1
                continue
            started = True
        if col > placed.last_col:
            break

        if col in placed.clip:
            kinds.append('clip')
            chars.append(placed.clip[col])
        elif col in placed.anchor:
            kinds.append('anchor')
            chars.append(placed.anchor[col])
        else:
            kinds.append('gap')
            chars.append(GAP)

    # Collapse consecutive columns of the same kind into one span.
    runs: List[Tuple[str, str]] = []
    for kind, char in zip(kinds, chars):
        if runs and runs[-1][0] == kind:
            runs[-1] = (kind, runs[-1][1] + char)
        else:
            runs.append((kind, char))

    return lead, runs


def render_reference(
    reference: str,
    reference_offsets: Sequence[int],
    columns: Sequence[Column],
) -> Tuple[str, int]:
    """
    Render the contig sequence over the shared columns.

    Parameters
    ----------
    reference : str
        Contig sequence covering ``reference_offsets``, in that order.
    reference_offsets : sequence of int
        Offsets the reference occupies.
    columns : sequence of Column
        Shared column order.

    Returns
    -------
    tuple
        ``(text, lead)``. Insertion columns show ``-``, since the reference has
        no base there by definition.
    """
    if not reference or not reference_offsets:
        return '', 0

    by_offset = dict(zip(reference_offsets, reference))
    first, last = min(reference_offsets), max(reference_offsets)

    lead = 0
    seen = False
    chars: List[str] = []

    for offset, ins in columns:
        if not seen:
            if offset < first:
                lead += 1
                continue
            seen = True
        if offset > last:
            break
        chars.append(by_offset.get(offset, GAP) if ins == 0 else GAP)

    return ''.join(chars), lead


def terminus_column(columns: Sequence[Column], is_left: bool) -> int:
    """
    Find the rendered column index of the contig/overhang boundary.

    Parameters
    ----------
    columns : sequence of Column
        Shared column order.
    is_left : bool
        True for a 5'-terminus overhang.

    Returns
    -------
    int
        Index of the boundary. At the left end this is the left edge of the
        terminal base; at the right end the edge after it, including any
        insertion columns that attach to it.
    """
    for i, (offset, ins) in enumerate(columns):
        if offset == 0 and ins == 0:
            if is_left:
                return i
            # At the right end the boundary sits after the terminal base and
            # after anything inserted at it.
            j = i + 1
            while j < len(columns) and columns[j][0] == 0:
                j += 1
            return j
    return 0
