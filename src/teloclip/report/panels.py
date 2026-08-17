"""
Per-read alignment panels for the ``teloclip extend`` HTML report.

This is the evidence half of the report: each overhang read drawn against the
contig terminus it supports, with the anchored portion and the soft clip
distinguished, so a reader can judge an extension rather than take it on
trust. Column geometry comes from :mod:`teloclip.report.layout`; this module
turns those columns into marked-up HTML.
"""

from dataclasses import dataclass
from html import escape
import re
from typing import Dict, List, Optional, Sequence, Tuple

from ..core.analysis import ContigStats, OverhangInfo
from .layout import (
    build_columns,
    place_read,
    render_reference,
    render_row,
    terminus_column,
)
from .text import fmt_delta_html, fmt_int


@dataclass(frozen=True)
class ReadRow:
    """One overhang read, laid out over the shared alignment columns."""

    read_name: str
    #: Blank columns before the read starts.
    lead: int
    #: ``(kind, text)`` runs, kind in {'clip', 'anchor', 'gap'}.
    runs: List[Tuple[str, str]]
    #: Bases the contig would gain from this read.
    net_gain: int
    #: Contig bases that would be trimmed to graft this clip on.
    trim: int
    #: True when this is the read the extension actually used.
    selected: bool
    #: True when the same read hangs off both ends of this contig.
    spans_both_ends: bool
    #: SAM record fields, for the hover table.
    cigar: str
    mapq: int
    flag: int
    aln_start: int
    aln_end: int
    clip_length: int
    anchor_length: int


@dataclass(frozen=True)
class EndPanel:
    """Everything needed to render one contig end."""

    contig: str
    end: str
    contig_length: int
    #: Reference row text and its leading blank columns.
    reference: str
    reference_lead: int
    rows: List[ReadRow]
    #: Rendered column index of the contig/overhang boundary.
    terminus_col: int
    #: ``(column_index, label)`` ticks for the position ruler.
    ticks: List[Tuple[int, str]]
    total_reads: int
    flagged: bool


def _anchor_window(
    overhangs: Sequence[OverhangInfo], contig_length: int, minimum: int, cap: int
) -> int:
    """
    Decide how much anchored reference to render at a contig end.

    At least ``minimum`` bases, or the longest anchor among the reads if that is
    longer, and never more than the contig itself.

    Parameters
    ----------
    overhangs : sequence of OverhangInfo
        Accepted overhangs at this end.
    contig_length : int
        Length of the contig.
    minimum : int
        Floor for the window.
    cap : int
        Ceiling, matching how much read sequence was retained.

    Returns
    -------
    int
        Window size in reference bases.
    """
    longest = max((oh.anchor_length for oh in overhangs), default=0)
    return max(1, min(contig_length, cap, max(minimum, longest)))


def _ruler(columns, terminus_col: int, step: int = 20) -> List[Tuple[int, str]]:
    """
    Build x-axis ticks with 0 at the contig terminus.

    Parameters
    ----------
    columns : sequence
        Shared column order, as ``(ref_offset, ins_index)`` pairs.
    terminus_col : int
        Rendered column index of the boundary.
    step : int, optional
        Spacing between ticks in reference bases (default: 20).

    Returns
    -------
    list of tuple
        ``(column_index, label)``, labelled by reference offset so insertion
        columns do not shift the scale.
    """
    ticks: List[Tuple[int, str]] = []
    for index, (offset, ins) in enumerate(columns):
        if ins != 0 or offset % step != 0:
            continue
        ticks.append((index, f'{offset:+d}' if offset else '0'))
    return ticks


def build_end_panel(
    contig: str,
    end: str,
    stats: ContigStats,
    overhangs: Sequence[OverhangInfo],
    terminal_sequence: str,
    selected_read: Optional[str],
    max_reads: int,
    max_overhang: int,
    flagged: bool,
    min_window: int = 500,
    window_cap: int = 1000,
) -> Optional[EndPanel]:
    """
    Lay out one contig end's overhang reads against the terminus.

    Each read is placed by walking its CIGAR, so indels shift it correctly
    rather than by a flat offset, and the column set is computed across all
    reads first so an insertion in one opens a gap in every other row.

    Parameters
    ----------
    contig : str
        Contig name.
    end : str
        Either ``'left'`` or ``'right'``.
    stats : ContigStats
        Statistics for the contig, for its length.
    overhangs : Sequence[OverhangInfo]
        Accepted overhangs at this end.
    terminal_sequence : str
        Terminal window of the original contig, oriented as in the assembly.
    selected_read : Optional[str]
        Name of the read the extension used, if any.
    max_reads : int
        Maximum reads to render, largest net gain first.
    max_overhang : int
        Maximum clipped bases to render per read.
    flagged : bool
        Whether this end was flagged for anomalous coverage.
    min_window : int, optional
        Minimum anchored reference bases to show (default: 500).
    window_cap : int, optional
        Maximum anchored reference bases to show (default: 1000).

    Returns
    -------
    Optional[EndPanel]
        The panel, or None when there are no overhangs to show.
    """
    if not overhangs:
        return None

    is_left = end == 'left'
    contig_length = stats.contig_length

    window = _anchor_window(overhangs, contig_length, min_window, window_cap)
    window = min(window, len(terminal_sequence) or window)

    ranked = sorted(overhangs, key=lambda oh: (oh.net_gain, oh.length), reverse=True)
    shown = ranked[:max_reads]

    placements = [
        place_read(
            cigar=oh.cigar,
            sequence=oh.read_seq,
            aln_start=oh.alignment_pos,
            aln_end=oh.alignment_end,
            contig_length=contig_length,
            clip_len=oh.clip_length,
            is_left=is_left,
            window=window,
            max_overhang=max_overhang,
            seq_offset=oh.read_seq_offset,
        )
        for oh in shown
    ]

    if is_left:
        ref_offsets = list(range(0, window))
        ref_seq = terminal_sequence[:window]
    else:
        ref_offsets = list(range(-window + 1, 1))
        ref_seq = terminal_sequence[-window:]

    columns = build_columns(placements, ref_offsets)
    ref_text, ref_lead = render_reference(ref_seq, ref_offsets, columns)

    rows: List[ReadRow] = []
    for oh, placed in zip(shown, placements):
        lead, runs = render_row(placed, columns)
        rows.append(
            ReadRow(
                read_name=oh.read_name,
                lead=lead,
                runs=runs,
                net_gain=oh.net_gain,
                trim=oh.length - oh.net_gain,
                selected=oh.read_name == selected_read,
                spans_both_ends=oh.spans_both_ends,
                cigar=oh.cigar,
                mapq=oh.mapq,
                flag=oh.flag,
                aln_start=oh.alignment_pos,
                aln_end=oh.alignment_end,
                clip_length=oh.clip_length,
                anchor_length=oh.anchor_length,
            )
        )

    return EndPanel(
        contig=contig,
        end=end,
        contig_length=contig_length,
        reference=ref_text,
        reference_lead=ref_lead,
        rows=rows,
        terminus_col=terminus_column(columns, is_left),
        ticks=_ruler(columns, terminus_column(columns, is_left)),
        total_reads=len(overhangs),
        flagged=flagged,
    )


def _render_sequence(
    seq: str, css_class: str, motif_pattern: Optional[re.Pattern]
) -> str:
    """
    Render a sequence as HTML, marking motif matches.

    Parameters
    ----------
    seq : str
        Sequence to render.
    css_class : str
        Class applied to the span.
    motif_pattern : Optional[re.Pattern]
        Compiled alternation of the motifs to highlight, or None.

    Returns
    -------
    str
        HTML markup.
    """
    if not seq:
        return ''
    if motif_pattern is None:
        return f'<span class="{css_class}">{escape(seq)}</span>'

    parts: List[str] = []
    cursor = 0
    for match in motif_pattern.finditer(seq):
        if match.start() > cursor:
            parts.append(escape(seq[cursor : match.start()]))
        parts.append(f'<mark>{escape(match.group(0))}</mark>')
        cursor = match.end()
    if cursor < len(seq):
        parts.append(escape(seq[cursor:]))
    return f'<span class="{css_class}">{"".join(parts)}</span>'


def _render_panel(panel: EndPanel, motif_pattern: Optional[re.Pattern]) -> str:
    """
    Render one contig end as a scrollable alignment block.

    Parameters
    ----------
    panel : EndPanel
        Prepared layout for this end.
    motif_pattern : Optional[re.Pattern]
        Motifs to highlight, or None.

    Returns
    -------
    str
        HTML markup.
    """
    is_left = panel.end == 'left'
    lines: List[str] = []

    # Position ruler, 0 at the contig terminus. Labels are anchored by their
    # left edge at the tick column and allowed to overflow, so they do not
    # widen the grid.
    ruler = ''.join(
        f'<span class="tick" style="left:{col}ch"><i>{escape(label)}</i></span>'
        for col, label in panel.ticks
    )
    lines.append(
        '<div class="aln-row aln-ruler">'
        '<span class="aln-label"></span>'
        f'<span class="aln-seq ruler">{ruler}</span>'
        '</div>'
    )

    lines.append(
        '<div class="aln-row aln-ref">'
        '<span class="aln-label" title="Original contig sequence">contig</span>'
        f'<span class="aln-seq" style="padding-left:{panel.reference_lead}ch">'
        f'{_render_sequence(panel.reference, "ref", motif_pattern)}</span>'
        '</div>'
    )

    for row in panel.rows:
        body = ''.join(
            _render_sequence(text, kind, motif_pattern if kind != 'gap' else None)
            for kind, text in row.runs
        )

        marks = ''
        if row.selected:
            marks += '<span class="pick" title="Used for this extension">&#9656;</span>'
        if row.spans_both_ends:
            marks += (
                '<span class="both" title="This read overhangs BOTH ends of '
                'the contig">&#8646;</span>'
            )

        label = row.read_name if len(row.read_name) <= 18 else row.read_name[:17] + '…'

        # Attributes for the hover table. Kept as data-* rather than a title so
        # the tooltip can be laid out as a table instead of one flat string.
        meta = (
            f' data-read="{escape(row.read_name)}"'
            f' data-flag="{row.flag}"'
            f' data-mapq="{"–" if row.mapq < 0 else row.mapq}"'
            f' data-cigar="{escape(row.cigar or "–")}"'
            f' data-span="{fmt_int(row.aln_start)}–{fmt_int(row.aln_end)}"'
            f' data-clip="{fmt_int(row.clip_length)}"'
            f' data-anchor="{fmt_int(row.anchor_length)}"'
            f' data-gain="{fmt_delta_html(row.net_gain)}"'
            f' data-trim="{fmt_int(row.trim)}"'
            f' data-end="{panel.end}"'
            f' data-both="{"yes" if row.spans_both_ends else "no"}"'
            f' data-selected="{"yes" if row.selected else "no"}"'
        )

        classes = 'aln-row aln-read'
        if row.selected:
            classes += ' is-selected'
        if row.spans_both_ends:
            classes += ' is-both'

        lines.append(
            f'<div class="{classes}"{meta}>'
            f'<span class="aln-label">{marks}{escape(label)}</span>'
            f'<span class="aln-seq" style="padding-left:{row.lead}ch">{body}</span>'
            '</div>'
        )

    hidden = panel.total_reads - len(panel.rows)
    note = (
        f'<p class="note">Showing the {len(panel.rows)} reads contributing most '
        f'sequence, of {panel.total_reads}. {hidden} not shown.</p>'
        if hidden > 0
        else ''
    )

    chips = f'<span class="chip chip-end">{panel.end} end</span>'
    if panel.flagged:
        chips += (
            '<span class="chip chip-flag" title="Overhang depth far above the '
            'rest of the assembly">&#9888; anomalous depth</span>'
        )
    if any(r.spans_both_ends for r in panel.rows):
        chips += (
            '<span class="chip chip-both" title="At least one read hangs off '
            'both ends of this contig">&#8646; spans both ends</span>'
        )

    return f"""
<section class="panel {'end-left' if is_left else 'end-right'}">
  <h4>{escape(panel.contig)} {chips}</h4>
  <div class="aln-scroll">
    <div class="aln" style="--terminus:{panel.terminus_col}ch">
      <div class="terminus" aria-hidden="true"></div>
      {''.join(lines)}
    </div>
  </div>
  {note}
</section>
"""


def render_contig_panels(
    contig: str,
    stats: ContigStats,
    terminal_sequences: Tuple[str, str],
    selected_reads: Dict[str, Optional[str]],
    flagged_left: bool,
    flagged_right: bool,
    motifs: Sequence[str] = (),
    max_reads: int = 25,
    max_overhang: int = 300,
    min_window: int = 500,
    window_cap: int = 1000,
) -> List[str]:
    """
    Render both alignment blocks for one contig.

    Exposed separately so the caller can render while a contig's read sequences
    are still in memory and drop them immediately after, rather than retaining
    them for the whole assembly.

    Parameters
    ----------
    contig : str
        Contig name.
    stats : ContigStats
        Overhang statistics for the contig.
    terminal_sequences : tuple of str
        ``(left_window, right_window)`` of the original contig.
    selected_reads : dict
        ``{'left': read, 'right': read}`` for applied extensions.
    flagged_left : bool
        Whether the 5' end was flagged for anomalous coverage.
    flagged_right : bool
        Whether the 3' end was flagged.
    motifs : sequence of str, optional
        Motifs to highlight.
    max_reads : int, optional
        Maximum reads rendered per end (default: 25).
    max_overhang : int, optional
        Maximum clipped bases rendered per read (default: 300).
    min_window : int, optional
        Minimum anchored reference bases shown (default: 500).
    window_cap : int, optional
        Maximum anchored reference bases shown (default: 1000).

    Returns
    -------
    list of str
        Rendered HTML blocks, one per end that had overhangs.
    """
    pattern = motif_pattern(motifs)
    left_window, right_window = terminal_sequences

    blocks: List[str] = []
    for end, overhangs, window, flagged in (
        ('left', stats.left_overhangs, left_window, flagged_left),
        ('right', stats.right_overhangs, right_window, flagged_right),
    ):
        panel = build_end_panel(
            contig=contig,
            end=end,
            stats=stats,
            overhangs=overhangs,
            terminal_sequence=window,
            selected_read=(selected_reads or {}).get(end),
            max_reads=max_reads,
            max_overhang=max_overhang,
            flagged=flagged,
            min_window=min_window,
            window_cap=window_cap,
        )
        if panel is not None:
            blocks.append(_render_panel(panel, pattern))
    return blocks


def motif_pattern(motifs: Sequence[str]) -> Optional[re.Pattern]:
    """
    Compile the motif alternation used for highlighting.

    Parameters
    ----------
    motifs : sequence of str
        Motif strings.

    Returns
    -------
    Optional[re.Pattern]
        Compiled pattern, or None when no motifs were requested.
    """
    if not motifs:
        return None
    return re.compile(
        '|'.join(re.escape(m) for m in sorted(set(motifs), key=len, reverse=True)),
        re.IGNORECASE,
    )
