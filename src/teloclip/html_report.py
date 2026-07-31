"""
Self-contained HTML report for ``teloclip extend``.

The Markdown report answers "what changed". This one answers "should I believe
it", by showing the evidence: every overhang read laid out against the contig
terminus it supports, and the distribution of overhang depth across the
assembly so anomalies are visible rather than merely listed.

The output is a single file with no external references — no CDN scripts, no
web fonts, no separate stylesheet — so it survives being emailed, copied onto a
cluster, or opened years later.
"""

from dataclasses import dataclass
from html import escape
import re
from typing import Dict, List, Optional, Sequence, Tuple

from .alignment_layout import (
    build_columns,
    place_read,
    render_reference,
    render_row,
    terminus_column,
)
from .analysis import ContigStats, OverhangInfo

# --- Palette -----------------------------------------------------------------
#
# Validated categorical slots 1 and 2 (blue, orange) carry the only identity
# encoding in the report: which contig end an overhang belongs to. Two series
# clear every gate on the all-pairs list in both modes.
#
# Flagged contigs are NOT given a third hue. Status is carried by a ring, a
# label and an icon, so identity and status never compete for the same channel.

_LIGHT = {
    'surface': '#fcfcfb',
    'plane': '#f9f9f7',
    'ink': '#0b0b0b',
    'ink2': '#52514e',
    'muted': '#898781',
    'grid': '#e1e0d9',
    'axis': '#c3c2b7',
    'left': '#2a78d6',
    'right': '#eb6834',
    'critical': '#d03b3b',
    'motif': 'rgba(12, 163, 12, 0.18)',
}

_DARK = {
    'surface': '#1a1a19',
    'plane': '#0d0d0d',
    'ink': '#ffffff',
    'ink2': '#c3c2b7',
    'muted': '#898781',
    'grid': '#2c2c2a',
    'axis': '#383835',
    'left': '#3987e5',
    'right': '#d95926',
    'critical': '#d03b3b',
    'motif': 'rgba(12, 163, 12, 0.30)',
}


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


def _fmt(value: int) -> str:
    """
    Format an integer with thousands separators.

    Parameters
    ----------
    value : int
        Number to format.

    Returns
    -------
    str
        Formatted number.
    """
    return f'{value:,}'


def _delta(value: int) -> str:
    """
    Format a signed change, using an en dash for zero.

    Parameters
    ----------
    value : int
        Change to format.

    Returns
    -------
    str
        Formatted change.
    """
    return '–' if value == 0 else f'{value:+,}'


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
            f' data-span="{_fmt(row.aln_start)}–{_fmt(row.aln_end)}"'
            f' data-clip="{_fmt(row.clip_length)}"'
            f' data-anchor="{_fmt(row.anchor_length)}"'
            f' data-gain="{_delta(row.net_gain)}"'
            f' data-trim="{_fmt(row.trim)}"'
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
    pattern = _motif_pattern(motifs)
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


def _motif_pattern(motifs: Sequence[str]) -> Optional[re.Pattern]:
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


def _coverage_chart(
    stats_dict: Dict[str, ContigStats],
    flagged_left: Sequence[str],
    flagged_right: Sequence[str],
) -> Tuple[str, str]:
    """
    Draw overhang depth per contig end as an inline SVG strip plot.

    The question is "does any end stand out", so the form is a strip plot with
    a median reference line: distribution first, individual values on hover.
    A bar chart would imply the contig order means something, and would not
    survive an assembly with hundreds of contigs.

    Parameters
    ----------
    stats_dict : Dict[str, ContigStats]
        Per-contig overhang statistics.
    flagged_left : Sequence[str]
        Contigs flagged at the left end.
    flagged_right : Sequence[str]
        Contigs flagged at the right end.

    Returns
    -------
    Tuple[str, str]
        The chart markup and the equivalent table markup.
    """
    names = list(stats_dict)
    if not names:
        return '', ''

    left_flags, right_flags = set(flagged_left), set(flagged_right)
    points = []
    for i, name in enumerate(names):
        s = stats_dict[name]
        points.append((i, name, 'left', s.left_count, name in left_flags))
        points.append((i, name, 'right', s.right_count, name in right_flags))

    counts = [p[3] for p in points]
    peak = max(counts) or 1
    ordered = sorted(counts)
    mid = len(ordered) // 2
    median = ordered[mid] if len(ordered) % 2 else (ordered[mid - 1] + ordered[mid]) / 2

    # Geometry. The plot fills the card at a comfortable width and only grows
    # (and scrolls) once there are enough contigs that points would collide.
    # The top pad has to clear the direct label sitting above the tallest
    # point, or a flagged contig at the peak has its name cropped.
    pad_l, pad_r, pad_t, pad_b = 52, 16, 30, 44
    target_w = 880
    step = max(16, target_w // max(1, len(names)))
    plot_w = max(target_w, step * max(1, len(names)))
    plot_h = 230
    width, height = plot_w + pad_l + pad_r, plot_h + pad_t + pad_b

    def x_of(i: int, end: str) -> float:
        """
        Horizontal position for a contig index, jittered by end.

        Parameters
        ----------
        i : int
            Contig index.
        end : str
            Either 'left' or 'right'.

        Returns
        -------
        float
            X coordinate.
        """
        nudge = -step * 0.16 if end == 'left' else step * 0.16
        return pad_l + step * (i + 0.5) + nudge

    def y_of(count: int) -> float:
        """
        Vertical position for a count.

        Parameters
        ----------
        count : int
            Overhang read count.

        Returns
        -------
        float
            Y coordinate.
        """
        return pad_t + plot_h - (count / peak) * plot_h

    svg: List[str] = []

    # Recessive solid hairline grid, four bands.
    for frac in (0, 0.25, 0.5, 0.75, 1.0):
        value = round(peak * frac)
        y = y_of(value)
        svg.append(
            f'<line class="grid" x1="{pad_l}" y1="{y:.1f}" '
            f'x2="{pad_l + plot_w}" y2="{y:.1f}"/>'
        )
        svg.append(
            f'<text class="tick" x="{pad_l - 8}" y="{y + 3.5:.1f}" '
            f'text-anchor="end">{_fmt(value)}</text>'
        )

    # Median reference: the line the flagged points deviate from.
    y_med = y_of(median)
    svg.append(
        f'<line class="median" x1="{pad_l}" y1="{y_med:.1f}" '
        f'x2="{pad_l + plot_w}" y2="{y_med:.1f}"/>'
    )
    svg.append(
        f'<text class="median-label" x="{pad_l + plot_w - 4}" '
        f'y="{y_med - 6:.1f}" text-anchor="end">median {median:g}</text>'
    )

    for i, name, end, count, is_flagged in points:
        cx, cy = x_of(i, end), y_of(count)
        classes = f'pt pt-{end}' + (' pt-flag' if is_flagged else '')
        # A generous transparent hit target sits under the visible mark.
        svg.append(
            f'<g class="{classes}" tabindex="0" role="listitem"'
            f' data-contig="{escape(name)}" data-end="{end}"'
            f' data-depth="{_fmt(count)}"'
            f' data-flagged="{"yes" if is_flagged else "no"}"'
            f' data-x="{cx:.1f}" data-y="{cy:.1f}">'
            f'<title>{escape(name)} — {end} end — {_fmt(count)} reads</title>'
            f'<circle class="hit" cx="{cx:.1f}" cy="{cy:.1f}" r="12"/>'
            f'<circle class="dot" cx="{cx:.1f}" cy="{cy:.1f}" r="4.5"/>'
            + (
                f'<circle class="ring" cx="{cx:.1f}" cy="{cy:.1f}" r="8.5"/>'
                if is_flagged
                else ''
            )
            + '</g>'
        )

    # Direct-label only the flagged points; everything else is on hover.
    for i, name, end, count, is_flagged in points:
        if not is_flagged:
            continue
        # Keep the label inside the plot: a flagged contig at either extreme
        # would otherwise have its name cropped by the card edge.
        half = 3.5 * len(name)
        lx = min(max(x_of(i, end), pad_l + half), pad_l + plot_w - half)
        svg.append(
            f'<text class="pt-label" x="{lx:.1f}" '
            f'y="{y_of(count) - 13:.1f}" text-anchor="middle">{escape(name)}</text>'
        )

    svg.append(
        f'<line class="axis" x1="{pad_l}" y1="{pad_t + plot_h}" '
        f'x2="{pad_l + plot_w}" y2="{pad_t + plot_h}"/>'
    )
    svg.append(
        f'<text class="axis-title" x="{pad_l}" y="{height - 8}">'
        f'{len(names)} contigs, in assembly order</text>'
    )

    chart = f"""
<div class="chart-scroll">
  <svg viewBox="0 0 {width} {height}" style="min-width:{width}px;width:100%"
       height="{height}" preserveAspectRatio="xMinYMid meet"
       role="img" aria-label="Overhang read depth for each contig end">
    {''.join(svg)}
  </svg>
</div>
"""

    rows = ''.join(
        f'<tr{" class=is-flagged" if f else ""}><td>{escape(n)}</td>'
        f'<td class="num">{_fmt(c)}</td><td>{e}</td>'
        f'<td>{"anomalous" if f else "–"}</td></tr>'
        for _, n, e, c, f in points
    )
    table = f"""
<details class="table-view">
  <summary>Table view</summary>
  <table>
    <thead><tr><th>Contig</th><th class="num">Overhang reads</th>
    <th>End</th><th>Depth</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
</details>
"""
    return chart, table


def _css() -> str:
    """
    Build the report stylesheet.

    Returns
    -------
    str
        CSS text, including the dark-mode variants.
    """

    def block(colors: Dict[str, str]) -> str:
        """
        Emit custom properties for one theme.

        Parameters
        ----------
        colors : Dict[str, str]
            Role to colour mapping.

        Returns
        -------
        str
            CSS declarations.
        """
        return '\n'.join(f'  --{k}: {v};' for k, v in colors.items())

    return f"""
:root {{
  color-scheme: light;
{block(_LIGHT)}
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    color-scheme: dark;
{block(_DARK)}
  }}
}}
:root[data-theme="dark"] {{
  color-scheme: dark;
{block(_DARK)}
}}

* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  padding: 2rem 1.25rem 4rem;
  background: var(--plane);
  color: var(--ink);
  font: 15px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
        "Helvetica Neue", Arial, sans-serif;
}}
main {{ max-width: 1180px; margin: 0 auto; }}
h1 {{ font-size: 1.6rem; margin: 0 0 .25rem; letter-spacing: -.01em; }}
h2 {{ font-size: 1.15rem; margin: 2.5rem 0 .75rem; letter-spacing: -.005em; }}
h4 {{ font-size: .95rem; margin: 0 0 .6rem; font-weight: 600; }}
p.sub {{ color: var(--ink2); margin: 0 0 1.5rem; }}
p.note {{ color: var(--muted); font-size: .82rem; margin: .5rem 0 0; }}
a {{ color: var(--left); }}

.tiles {{
  display: grid; gap: .75rem;
  grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
}}
.tile {{
  background: var(--surface); border: 1px solid var(--grid);
  border-radius: 10px; padding: .85rem 1rem;
}}
.tile .k {{ color: var(--ink2); font-size: .78rem; text-transform: uppercase;
  letter-spacing: .04em; }}
.tile .v {{ font-size: 1.5rem; font-weight: 600; margin-top: .15rem; }}

table {{ border-collapse: collapse; width: 100%; font-size: .88rem; }}
th, td {{ text-align: left; padding: .4rem .6rem;
  border-bottom: 1px solid var(--grid); }}
th {{ color: var(--ink2); font-weight: 600; }}
td.num, th.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
tr.is-flagged td {{ color: var(--critical); }}
.scroll-x {{ overflow-x: auto; }}

.chart-scroll, .aln-scroll {{ overflow-x: auto; }}

.card {{
  background: var(--surface); border: 1px solid var(--grid);
  border-radius: 10px; padding: 1rem;
}}
svg .grid {{ stroke: var(--grid); stroke-width: 1; }}
svg .axis {{ stroke: var(--axis); stroke-width: 1; }}
svg .median {{ stroke: var(--axis); stroke-width: 2; }}
svg .tick, svg .axis-title, svg .median-label {{
  fill: var(--muted); font-size: 11px; font-variant-numeric: tabular-nums;
}}
svg .pt-label {{ fill: var(--ink2); font-size: 11px; font-weight: 600; }}
svg .hit {{ fill: transparent; }}
svg .dot {{ stroke: var(--surface); stroke-width: 2; }}
svg .pt-left .dot {{ fill: var(--left); }}
svg .pt-right .dot {{ fill: var(--right); }}
svg .ring {{ fill: none; stroke: var(--critical); stroke-width: 2; }}
svg .pt:hover .dot {{ r: 6; }}

.legend {{ display: flex; gap: 1.1rem; flex-wrap: wrap;
  margin: .75rem 0 0; font-size: .82rem; color: var(--ink2); }}
.legend span.sw {{ display: inline-block; width: 10px; height: 10px;
  border-radius: 50%; margin-right: .35rem; vertical-align: -1px; }}
.sw-left {{ background: var(--left); }}
.sw-right {{ background: var(--right); }}
.sw-flag {{ background: transparent; border: 2px solid var(--critical); }}

.panel {{
  background: var(--surface); border: 1px solid var(--grid);
  border-radius: 10px; padding: 1rem; margin: 0 0 .85rem;
}}
.chip {{
  display: inline-block; font-size: .72rem; font-weight: 600;
  padding: .1rem .45rem; border-radius: 999px; vertical-align: 2px;
}}
.chip-end {{ background: var(--grid); color: var(--ink2); }}
.end-left .chip-end {{ color: var(--left); }}
.end-right .chip-end {{ color: var(--right); }}
.chip-flag {{ color: var(--critical); border: 1px solid var(--critical); }}

/* The monospace font lives on the container so that the `ch` unit resolves
   identically for the sequence padding and for the terminus rule. Setting it
   only on the rows would leave the rule positioned in the body font's `ch`,
   which is a different width, and the marker would drift from the sequence. */
.aln {{
  position: relative; width: max-content; min-width: 100%;
  font: 12px/1.45 ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  /* Sticky read-name gutter, sized in the same `ch` grid as the sequences. */
  --gutter: 22ch;
}}
.terminus {{
  position: absolute; top: 0; bottom: 0;
  left: calc(var(--gutter) + var(--terminus)); width: 2px;
  background: var(--critical); opacity: .55;
}}
.aln-row {{ white-space: pre; }}
.aln-row.is-selected {{ background: color-mix(in srgb, var(--grid) 45%, transparent); }}
.aln-label {{
  display: inline-block; width: var(--gutter); padding-right: .6rem;
  font-weight: 600; color: var(--ink2); text-align: right;
  position: sticky; left: 0; background: var(--surface);
  overflow: hidden; text-overflow: ellipsis; white-space: nowrap;
}}
.aln-ref .aln-label {{ color: var(--ink); }}
.ref {{ color: var(--ink); }}
.anchor {{ color: var(--muted); }}
.gap {{ color: var(--axis); }}
.aln-ruler {{ position: relative; height: 1.45em; }}
.aln-ruler .aln-seq {{ position: relative; display: inline-block; }}
.ruler .tick {{
  /* No font-size here: `left` is expressed in `ch`, which must resolve
     against the 12px monospace grid the sequences use. Shrinking the tick
     would shrink its `ch` and drift the scale off the bases. */
  position: absolute; top: 0; white-space: nowrap;
  border-left: 1px solid var(--axis); padding-left: 2px;
}}
.ruler .tick i {{ font-style: normal; font-size: 10px; color: var(--muted); }}
.aln-read {{ cursor: default; }}
.aln-read:hover {{ background: color-mix(in srgb, var(--left) 10%, transparent); }}
.both {{ color: var(--critical); margin-right: .25rem; }}
.chip-both {{ color: var(--critical); border: 1px solid var(--critical); }}
.is-both .aln-label {{ color: var(--critical); }}

#tip {{
  position: fixed; z-index: 50; pointer-events: none; opacity: 0;
  transition: opacity .08s; background: var(--surface); color: var(--ink);
  border: 1px solid var(--axis); border-radius: 8px; padding: .5rem .6rem;
  box-shadow: 0 6px 24px rgba(0,0,0,.18); font-size: .78rem; max-width: 32rem;
}}
#tip.show {{ opacity: 1; }}
#tip table {{ border-collapse: collapse; width: auto; font-size: .78rem; }}
#tip td {{ border: 0; padding: .1rem .5rem .1rem 0; vertical-align: top; }}
#tip td.k {{ color: var(--ink2); white-space: nowrap; }}
#tip td.v {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  word-break: break-all; }}
#tip .tip-head {{ font-weight: 600; margin-bottom: .3rem; word-break: break-all; }}
svg .pt {{ cursor: pointer; }}
svg .pt-pin {{ fill: var(--ink2); font-size: 11px; font-weight: 600; }}
.end-left .clip {{ color: var(--left); font-weight: 600; }}
.end-right .clip {{ color: var(--right); font-weight: 600; }}
mark {{ background: var(--motif); color: inherit; border-radius: 2px; }}
.pick {{ color: var(--ink); margin-right: .25rem; }}

.table-view {{ margin-top: .75rem; font-size: .85rem; }}
.table-view summary {{ cursor: pointer; color: var(--ink2); }}
.table-view table {{ margin-top: .5rem; }}
footer {{ margin-top: 3rem; color: var(--muted); font-size: .8rem; }}
"""


def render_html_report(
    *,
    stats_dict: Dict[str, ContigStats],
    extensions_applied: Dict[str, dict],
    anomalous: Dict[str, List[str]],
    excluded_contigs: Sequence[str],
    warnings: Sequence[str],
    panels: Sequence[str],
    total_contigs: int,
    dry_run: bool,
    motifs: Sequence[str] = (),
    max_reads: int = 25,
    max_overhang: int = 300,
    version: str = '',
    command: str = '',
) -> str:
    """
    Render the full HTML report.

    Parameters
    ----------
    stats_dict : Dict[str, ContigStats]
        Overhang statistics for contigs with support.
    extensions_applied : Dict[str, dict]
        Extension records keyed by contig.
    anomalous : Dict[str, List[str]]
        Contigs flagged for anomalous coverage, keyed ``left_outliers`` and
        ``right_outliers``.
    excluded_contigs : Sequence[str]
        Contigs the user excluded.
    warnings : Sequence[str]
        Accumulated warnings.
    panels : Sequence[str]
        Alignment blocks already rendered by :func:`render_contig_panels`.
    total_contigs : int
        Contigs in the assembly.
    dry_run : bool
        Whether extensions were simulated.
    motifs : Sequence[str], optional
        Motif strings to highlight in sequence blocks.
    max_reads : int, optional
        Maximum reads rendered per contig end (default: 25).
    max_overhang : int, optional
        Maximum clipped bases rendered per read (default: 300).
    version : str, optional
        Teloclip version, for the footer.
    command : str, optional
        Command line that produced the report.

    Returns
    -------
    str
        A complete, self-contained HTML document.
    """
    motif_pattern = _motif_pattern(motifs)

    # --- Summary -------------------------------------------------------------
    ends_extended = 0
    net_gained = 0
    trimmed = 0
    for info in extensions_applied.values():
        for side in ('left', 'right'):
            if not info.get(f'has_{side}_extension'):
                continue
            ends_extended += 1
            gain = info.get(f'{side}_net_gain')
            if gain is None:
                gain = info.get(f'{side}_overhang_length', 0) - info.get(
                    f'{side}_trim_length', 0
                )
            net_gained += gain
            trimmed += info.get(f'{side}_trim_length', 0)

    flagged_all = sorted(
        set(anomalous.get('left_outliers', []))
        | set(anomalous.get('right_outliers', []))
    )

    tiles = [
        ('Mode', 'dry run' if dry_run else 'applied'),
        ('Contigs', _fmt(total_contigs)),
        ('Contigs extended', _fmt(len(extensions_applied))),
        ('Ends extended', f'{_fmt(ends_extended)} of {_fmt(total_contigs * 2)}'),
        ('Net bases gained', _fmt(net_gained)),
        ('Bases trimmed', _fmt(trimmed)),
        ('Flagged contigs', _fmt(len(flagged_all))),
        ('Warnings', _fmt(len(warnings))),
    ]
    tiles_html = ''.join(
        f'<div class="tile"><div class="k">{escape(k)}</div>'
        f'<div class="v">{escape(str(v))}</div></div>'
        for k, v in tiles
    )

    # --- Extensions table ----------------------------------------------------
    ext_rows = []
    for contig in sorted(extensions_applied):
        info = extensions_applied[contig]

        def net(side: str, info=info) -> int:
            """
            Net gain at one end of this contig.

            Parameters
            ----------
            side : str
                Either 'left' or 'right'.
            info : dict
                Extension record.

            Returns
            -------
            int
                Bases gained.
            """
            if not info.get(f'has_{side}_extension'):
                return 0
            value = info.get(f'{side}_net_gain')
            if value is None:
                value = info.get(f'{side}_overhang_length', 0) - info.get(
                    f'{side}_trim_length', 0
                )
            return value

        left, right = net('left'), net('right')
        ext_rows.append(
            f'<tr><td>{escape(contig)}</td>'
            f'<td class="num">{_fmt(info.get("original_length", 0))}</td>'
            f'<td class="num">{_fmt(info.get("final_length", 0))}</td>'
            f'<td class="num">{_delta(left)}</td>'
            f'<td class="num">{_delta(right)}</td>'
            f'<td class="num">{_delta(left + right)}</td>'
            f'<td class="num">{_fmt(info.get("left_trim_length", 0))}</td>'
            f'<td class="num">{_fmt(info.get("right_trim_length", 0))}</td></tr>'
        )

    ext_table = (
        f"""
<div class="card scroll-x">
<table>
  <thead><tr><th>Contig</th><th class="num">Original bp</th>
  <th class="num">Final bp</th><th class="num">Left +bp</th>
  <th class="num">Right +bp</th><th class="num">Total +bp</th>
  <th class="num">Left trim</th><th class="num">Right trim</th></tr></thead>
  <tbody>{''.join(ext_rows)}</tbody>
</table>
</div>
<p class="note">The <code>+bp</code> columns are net of trimming, so
<code>Original bp + Total +bp = Final bp</code> on every row.</p>
"""
        if ext_rows
        else '<p class="note">No contigs were extended.</p>'
    )

    # --- Coverage chart ------------------------------------------------------
    chart, chart_table = _coverage_chart(
        stats_dict,
        anomalous.get('left_outliers', []),
        anomalous.get('right_outliers', []),
    )

    # --- Alignment panels ----------------------------------------------------
    panels_html = ''.join(panels) or (
        '<p class="note">No overhang reads were retained.</p>'
    )

    # --- Secondary sections --------------------------------------------------
    anomaly_section = ''
    if flagged_all:
        items = ''.join(f'<li><code>{escape(c)}</code></li>' for c in flagged_all)
        anomaly_section = f"""
<h2>Contigs with anomalous overhang coverage</h2>
<div class="card">
<p>These ends carry far more clipped reads than the rest of the assembly, which
usually indicates a collapsed repeat, an rDNA array, or an organellar contig
attracting reads from elsewhere. Extending from a single read is rarely
meaningful in those cases.</p>
<p><strong>They have not been excluded.</strong> Review them and re-run with
<code>--exclude-contigs</code> if you agree.</p>
<ul>{items}</ul>
</div>
"""

    excluded_section = ''
    if excluded_contigs:
        items = ''.join(
            f'<li><code>{escape(c)}</code></li>' for c in sorted(excluded_contigs)
        )
        excluded_section = f"""
<h2>Excluded contigs</h2>
<div class="card"><p>Left untouched at your request, and written to the output
unmodified.</p><ul>{items}</ul></div>
"""

    warnings_section = ''
    if warnings:
        items = ''.join(f'<li>{escape(w)}</li>' for w in warnings)
        warnings_section = f"""
<h2>Warnings</h2>
<div class="card"><ul>{items}</ul></div>
"""

    legend = """
<div class="legend">
  <span><span class="sw sw-left"></span>left end</span>
  <span><span class="sw sw-right"></span>right end</span>
  <span><span class="sw sw-flag"></span>&#9888; anomalous depth</span>
</div>
"""

    footer_bits = [b for b in (f'teloclip {version}' if version else '', command) if b]
    footer = '<br>'.join(escape(b) for b in footer_bits)

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Teloclip extend report</title>
<style>{_css()}</style>
</head>
<body>
<main>
  <h1>Teloclip extend report</h1>
  <p class="sub">{
        'Dry run — no sequences were written.' if dry_run else 'Extensions applied.'
    }</p>

  <h2>Summary</h2>
  <div class="tiles">{tiles_html}</div>

  <h2>Extensions</h2>
  {ext_table}

  <h2>Overhang depth across contigs</h2>
  <div class="card">
    {chart}
    {legend}
    <p class="note">Each point is one contig end. Hover for the contig name and
    depth. Points far above the median are worth investigating.</p>
    {chart_table}
  </div>

  <h2>Overhang alignments</h2>
  <p class="note">Reads are laid out against the contig terminus, marked by the
  vertical rule. Grey is the anchored portion of the read; colour is the soft
  clip. Hover a row for the read name and what it contributes.
  {'Highlighted bases are motif matches.' if motif_pattern else ''}</p>
  {panels_html}

  {anomaly_section}
  {excluded_section}
  {warnings_section}

  <footer>{footer}</footer>
</main>
<div id="tip" role="tooltip" aria-hidden="true"></div>
<script>
// Everything below is progressive enhancement. Without JavaScript the report
// still renders: the blocks simply start at their left edge, and every row and
// data point keeps a native title attribute.
(function () {{
  // Open each alignment on the contig/overhang junction rather than at the far
  // end of the longest clip, which is what a left-aligned scroll box shows.
  for (const box of document.querySelectorAll('.aln-scroll')) {{
    const rule = box.querySelector('.terminus');
    if (!rule) continue;
    const label = box.querySelector('.aln-label');
    const gutter = label ? label.getBoundingClientRect().width : 0;
    box.scrollLeft = rule.offsetLeft - gutter - (box.clientWidth - gutter) * 0.55;
  }}

  const tip = document.getElementById('tip');
  const esc = (v) => String(v).replace(/[&<>"]/g, (c) => (
    {{'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}}[c]));
  const table = (rows) => '<table>' + rows.map(
    ([k, v]) => '<tr><td class="k">' + esc(k) + '</td><td class="v">' +
                esc(v) + '</td></tr>').join('') + '</table>';

  function show(html, x, y) {{
    tip.innerHTML = html;
    tip.classList.add('show');
    tip.setAttribute('aria-hidden', 'false');
    const r = tip.getBoundingClientRect();
    // Keep the tooltip on screen near the cursor.
    const left = Math.min(Math.max(8, x + 14), window.innerWidth - r.width - 8);
    const top = y + r.height + 24 > window.innerHeight ? y - r.height - 12 : y + 18;
    tip.style.left = left + 'px';
    tip.style.top = Math.max(8, top) + 'px';
  }}
  function hide() {{
    tip.classList.remove('show');
    tip.setAttribute('aria-hidden', 'true');
  }}

  // Read rows: a table of the SAM record and what the read contributes.
  for (const row of document.querySelectorAll('.aln-read')) {{
    const d = row.dataset;
    const rows = [
      ['End', d.end],
      ['Alignment', d.span],
      ['CIGAR', d.cigar],
      ['FLAG', d.flag],
      ['MAPQ', d.mapq],
      ['Anchor', d.anchor + ' bp'],
      ['Soft clip', d.clip + ' bp'],
      ['Contig gains', d.gain + ' bp'],
      ['Contig trimmed', d.trim + ' bp'],
    ];
    if (d.both === 'yes') rows.push(['Note', 'overhangs both ends of this contig']);
    if (d.selected === 'yes') rows.push(['Note', 'used for this extension']);
    const html = '<div class="tip-head">' + esc(d.read) + '</div>' + table(rows);
    row.addEventListener('mousemove', (e) => show(html, e.clientX, e.clientY));
    row.addEventListener('mouseleave', hide);
  }}

  // Chart points: contig name on hover, and on click, pinned as a label so it
  // survives moving the mouse away.
  const svg = document.querySelector('.chart-scroll svg');
  for (const pt of document.querySelectorAll('svg g.pt')) {{
    const d = pt.dataset;
    const rows = [['Contig', d.contig], ['End', d.end],
                  ['Overhang reads', d.depth]];
    if (d.flagged === 'yes') rows.push(['Depth', 'anomalous']);
    const html = table(rows);
    pt.addEventListener('mousemove', (e) => show(html, e.clientX, e.clientY));
    pt.addEventListener('mouseleave', hide);
    pt.addEventListener('focus', () => {{
      const r = pt.getBoundingClientRect();
      show(html, r.right, r.top);
    }});
    pt.addEventListener('blur', hide);
    pt.addEventListener('click', () => {{
      const existing = pt.querySelector('.pt-pin');
      if (existing) {{ existing.remove(); return; }}
      const t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      t.setAttribute('class', 'pt-pin');
      t.setAttribute('x', d.x);
      t.setAttribute('y', String(parseFloat(d.y) - 13));
      t.setAttribute('text-anchor', 'middle');
      t.textContent = d.contig;
      pt.appendChild(t);
    }});
  }}
  if (svg) svg.addEventListener('mouseleave', hide);
  window.addEventListener('scroll', hide, {{passive: true}});
}})();
</script>
</body>
</html>
"""
