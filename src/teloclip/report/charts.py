"""
Inline SVG charts for the ``teloclip extend`` HTML report.

Every chart here is emitted as hand-written SVG rather than drawn by a
plotting library, so a report stays a single self-contained file with nothing
to fetch. Each chart function returns a ``(chart, table)`` pair: the second
element is an equivalent ``<details>`` table, so the same information is
reachable without SVG, without colour, and by a screen reader.

Styling lives in :mod:`teloclip.report.css`; these functions emit class names
only.
"""

from html import escape
from typing import Dict, List, Sequence, Tuple

from ..core.analysis import ContigStats
from .text import fmt_int


def coverage_chart(
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
            f'text-anchor="end">{fmt_int(value)}</text>'
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
            f' data-depth="{fmt_int(count)}"'
            f' data-flagged="{"yes" if is_flagged else "no"}"'
            f' data-x="{cx:.1f}" data-y="{cy:.1f}">'
            f'<title>{escape(name)} — {end} end — {fmt_int(count)} reads</title>'
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
        f'<td class="num">{fmt_int(c)}</td><td>{e}</td>'
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
