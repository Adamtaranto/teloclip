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
import statistics
from typing import Dict, List, Sequence, Tuple

from ..core.analysis import ContigStats
from .text import fmt_int

# Shared plot geometry. Defined once so the three charts line up with each
# other down the page: a reader comparing them should not have to re-locate the
# plotting area each time.
PAD_L, PAD_R, PAD_T, PAD_B = 52, 16, 30, 44
TARGET_W = 880
PLOT_H = 230

# Below this many reads at an end, a density curve is a fiction drawn through
# too few points. Those ends are shown as individual marks instead.
MIN_READS_FOR_DENSITY = 4

# Bins across the shared length scale. Enough to show bimodality, few enough
# that a typical end with a few dozen reads does not turn into a comb.
DENSITY_BINS = 14

# The length axis is clipped here so that one read running far past the contig
# end cannot flatten every violin in the assembly into a line.
LENGTH_CLIP_QUANTILE = 0.99


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


def _quantile(sorted_values: List[float], q: float) -> float:
    """
    Linear-interpolated quantile of an already-sorted sequence.

    Written out rather than taken from ``statistics.quantiles`` so that the
    behaviour on very small samples is explicit: with one value that value is
    returned, and no exception is raised for n < 2.

    Parameters
    ----------
    sorted_values : List[float]
        Values in ascending order. Must not be empty.
    q : float
        Quantile in the range 0 to 1.

    Returns
    -------
    float
        The interpolated quantile.
    """
    if len(sorted_values) == 1:
        return float(sorted_values[0])

    position = q * (len(sorted_values) - 1)
    low = int(position)
    high = min(low + 1, len(sorted_values) - 1)
    weight = position - low
    return sorted_values[low] * (1 - weight) + sorted_values[high] * weight


def _density_profile(values: List[float], lo: float, hi: float) -> List[float]:
    """
    Bin values into a smoothed density profile over a fixed range.

    A histogram with a three-point moving average rather than a kernel density
    estimate. The smoothing is what stops the profile reading as a bar chart
    on its side, and avoids a dependency on scipy for a shape that is only ever
    read qualitatively.

    Parameters
    ----------
    values : List[float]
        Observations to bin.
    lo : float
        Lower edge of the shared scale.
    hi : float
        Upper edge of the shared scale.

    Returns
    -------
    List[float]
        One weight per bin, scaled so the largest is 1.0. All zeros when there
        is nothing to bin.
    """
    counts = [0.0] * DENSITY_BINS
    span = hi - lo

    for value in values:
        if span <= 0:
            index = 0
        else:
            index = int((value - lo) / span * DENSITY_BINS)
            index = min(max(index, 0), DENSITY_BINS - 1)
        counts[index] += 1

    # Two passes. One leaves a visibly ragged outline at the read counts these
    # ends typically carry, which reads as structure that is not in the data.
    smoothed = counts
    for _ in range(2):
        pass_result = []
        for i in range(DENSITY_BINS):
            window = smoothed[max(0, i - 1) : min(DENSITY_BINS, i + 2)]
            pass_result.append(sum(window) / len(window))
        smoothed = pass_result

    peak = max(smoothed)
    if peak == 0:
        return smoothed
    return [value / peak for value in smoothed]


def overhang_length_chart(
    stats_dict: Dict[str, ContigStats],
    flagged_left: Sequence[str],
    flagged_right: Sequence[str],
) -> Tuple[str, str]:
    """
    Draw the distribution of overhang lengths per contig end as split violins.

    Depth says how many reads support an end; this says how far they reach past
    it, which is what determines how much sequence an extension can recover. A
    tight distribution at a few hundred bases reads very differently from a
    broad one running into kilobases, and neither is visible in a count.

    The violin is split rather than paired: the left half of each shape is the
    contig's left end and the right half its right end, sharing one axis
    position. That keeps the two ends of a contig adjacent for comparison,
    which is the comparison a reader actually makes, and halves the width so an
    assembly with many contigs still fits.

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
        The chart markup and the equivalent table markup. Both are empty when
        there is nothing with an overhang to draw.
    """
    names = [
        name for name, s in stats_dict.items() if s.left_overhangs or s.right_overhangs
    ]
    if not names:
        return '', ''

    left_flags, right_flags = set(flagged_left), set(flagged_right)

    # One shared length scale, so violin heights are comparable between
    # contigs. Clipped at a high quantile because a single read running far
    # past the contig end would otherwise compress every other shape to a line.
    all_lengths = sorted(
        float(oh.length)
        for s in stats_dict.values()
        for oh in (*s.left_overhangs, *s.right_overhangs)
    )
    if not all_lengths:
        return '', ''

    hi = _quantile(all_lengths, LENGTH_CLIP_QUANTILE)
    longest = all_lengths[-1]
    # Only disclose clipping when the axis maximum a reader sees actually
    # differs from the longest overhang. On a tight distribution the 99th
    # percentile rounds to the same number as the maximum, and announcing a
    # clip there would put a caveat on every chart that did not need one.
    clipped = round(hi) < round(longest)
    if hi <= 0:
        hi = longest or 1.0

    step = max(24, TARGET_W // max(1, len(names)))
    plot_w = max(TARGET_W, step * max(1, len(names)))
    width = plot_w + PAD_L + PAD_R
    height = PLOT_H + PAD_T + PAD_B

    # Half-width of a violin at full density, leaving a gap between neighbours.
    half_w = step * 0.38

    def x_of(i: int) -> float:
        """
        Centre line for a contig's violin.

        Parameters
        ----------
        i : int
            Contig index.

        Returns
        -------
        float
            X coordinate of the shared centre line.
        """
        return PAD_L + step * (i + 0.5)

    def y_of(length: float) -> float:
        """
        Vertical position for an overhang length.

        Parameters
        ----------
        length : float
            Overhang length in bases.

        Returns
        -------
        float
            Y coordinate, clamped to the plotting area.
        """
        fraction = min(length / hi, 1.0) if hi else 0.0
        return PAD_T + PLOT_H - fraction * PLOT_H

    svg: List[str] = []

    for frac in (0, 0.25, 0.5, 0.75, 1.0):
        value = hi * frac
        y = y_of(value)
        svg.append(
            f'<line class="grid" x1="{PAD_L}" y1="{y:.1f}" '
            f'x2="{PAD_L + plot_w}" y2="{y:.1f}"/>'
        )
        svg.append(
            f'<text class="tick" x="{PAD_L - 8}" y="{y + 3.5:.1f}" '
            f'text-anchor="end">{fmt_int(round(value))}</text>'
        )

    for i, name in enumerate(names):
        stats = stats_dict[name]
        centre = x_of(i)

        for end, overhangs, flagged in (
            ('left', stats.left_overhangs, name in left_flags),
            ('right', stats.right_overhangs, name in right_flags),
        ):
            if not overhangs:
                continue

            lengths = sorted(float(oh.length) for oh in overhangs)
            # -1 draws the shape to the left of the centre line, +1 to the right.
            direction = -1 if end == 'left' else 1
            median = statistics.median(lengths)
            classes = f'vio vio-{end}' + (' vio-flag' if flagged else '')

            summary = (
                f'{escape(name)} — {end} end — {len(lengths)} reads, '
                f'median {fmt_int(round(median))} bp, '
                f'longest {fmt_int(round(lengths[-1]))} bp'
            )
            group = [
                f'<g class="{classes}" tabindex="0" role="listitem"'
                f' data-contig="{escape(name)}" data-end="{end}"'
                f' data-depth="{fmt_int(len(lengths))}"'
                f' data-median="{fmt_int(round(median))}"'
                f' data-longest="{fmt_int(round(lengths[-1]))}"'
                f' data-flagged="{"yes" if flagged else "no"}"'
                f'><title>{summary}</title>'
            ]

            if len(lengths) < MIN_READS_FOR_DENSITY:
                # A density curve through three points describes the estimator,
                # not the data. Show the reads themselves instead.
                for offset, length in enumerate(lengths):
                    jitter = direction * (half_w * 0.35 + (offset % 2) * 3.0)
                    group.append(
                        f'<circle class="vio-dot" cx="{centre + jitter:.1f}" '
                        f'cy="{y_of(length):.1f}" r="3.5"/>'
                    )
            else:
                profile = _density_profile(lengths, 0.0, hi)
                band = PLOT_H / DENSITY_BINS

                # Out along the density, then back down the centre line, so the
                # polygon closes flush against its opposite half. Closed at the
                # first and last *occupied* bin rather than at the plot edges:
                # anchoring to the edges stretched every shape over the full
                # height regardless of the range its reads actually covered.
                occupied = [b for b, weight in enumerate(profile) if weight > 0]
                first_bin, last_bin = occupied[0], occupied[-1]

                def y_of_bin(b: int, band=band) -> float:
                    """
                    Centre line of a density bin.

                    Parameters
                    ----------
                    b : int
                        Bin index, counting up from the bottom of the scale.
                    band : float
                        Height of one bin.

                    Returns
                    -------
                    float
                        Y coordinate of the bin's centre.
                    """
                    return PAD_T + PLOT_H - (b + 0.5) * band

                outer = [
                    f'{centre + direction * profile[b] * half_w:.1f},{y_of_bin(b):.1f}'
                    for b in range(first_bin, last_bin + 1)
                ]
                points = (
                    [f'{centre:.1f},{y_of_bin(first_bin):.1f}']
                    + outer
                    + [f'{centre:.1f},{y_of_bin(last_bin):.1f}']
                )
                group.append(f'<polygon class="vio-body" points="{" ".join(points)}"/>')

                # Interquartile range and median, drawn over the density as a
                # thin rule: the shape gives the distribution, this gives the
                # numbers a reader compares between contigs.
                q1, q3 = _quantile(lengths, 0.25), _quantile(lengths, 0.75)
                box_x = centre + direction * 2.0
                group.append(
                    f'<line class="vio-iqr" x1="{box_x:.1f}" y1="{y_of(q1):.1f}" '
                    f'x2="{box_x:.1f}" y2="{y_of(q3):.1f}"/>'
                )

            # A short tick, not a rule spanning the whole half. At full width it
            # was the heaviest mark in the chart, and on an end drawn as a few
            # dots it read as a bar of its own rather than as their summary.
            y_med = y_of(median)
            group.append(
                f'<line class="vio-median" x1="{centre:.1f}" y1="{y_med:.1f}" '
                f'x2="{centre + direction * half_w * 0.55:.1f}" y2="{y_med:.1f}"/>'
            )
            # Hit area spanning the full height of this half, so the shape does
            # not have to be landed on precisely.
            group.append(
                f'<rect class="hit" x="{min(centre, centre + direction * half_w):.1f}" '
                f'y="{PAD_T}" width="{half_w:.1f}" height="{PLOT_H}"/>'
            )
            group.append('</g>')
            svg.append(''.join(group))

    # Direct-label only the flagged ends, as the depth chart does.
    for i, name in enumerate(names):
        if name not in left_flags and name not in right_flags:
            continue
        half_label = 3.5 * len(name)
        lx = min(max(x_of(i), PAD_L + half_label), PAD_L + plot_w - half_label)
        svg.append(
            f'<text class="pt-label" x="{lx:.1f}" y="{PAD_T - 8}" '
            f'text-anchor="middle">{escape(name)}</text>'
        )

    svg.append(
        f'<line class="axis" x1="{PAD_L}" y1="{PAD_T + PLOT_H}" '
        f'x2="{PAD_L + plot_w}" y2="{PAD_T + PLOT_H}"/>'
    )
    axis_note = (
        f'{len(names)} contigs with overhangs, in assembly order. Overhang length in bp'
    )
    if clipped:
        axis_note += (
            f' (scale clipped at {fmt_int(round(hi))}; '
            f'longest is {fmt_int(round(longest))})'
        )
    svg.append(
        f'<text class="axis-title" x="{PAD_L}" y="{height - 8}">{axis_note}</text>'
    )

    chart = f"""
<div class="chart-scroll">
  <svg viewBox="0 0 {width} {height}" style="min-width:{width}px;width:100%"
       height="{height}" preserveAspectRatio="xMinYMid meet"
       role="img" aria-label="Distribution of overhang lengths for each contig end">
    {''.join(svg)}
  </svg>
</div>
"""

    rows = []
    for name in names:
        stats = stats_dict[name]
        for end, overhangs, flagged in (
            ('left', stats.left_overhangs, name in left_flags),
            ('right', stats.right_overhangs, name in right_flags),
        ):
            if not overhangs:
                continue
            lengths = sorted(float(oh.length) for oh in overhangs)
            rows.append(
                f'<tr{" class=is-flagged" if flagged else ""}>'
                f'<td>{escape(name)}</td><td>{end}</td>'
                f'<td class="num">{fmt_int(len(lengths))}</td>'
                f'<td class="num">{fmt_int(round(_quantile(lengths, 0.25)))}</td>'
                f'<td class="num">{fmt_int(round(statistics.median(lengths)))}</td>'
                f'<td class="num">{fmt_int(round(_quantile(lengths, 0.75)))}</td>'
                f'<td class="num">{fmt_int(round(lengths[-1]))}</td></tr>'
            )

    table = f"""
<details class="table-view">
  <summary>Table view</summary>
  <table>
    <thead><tr><th>Contig</th><th>End</th><th class="num">Reads</th>
    <th class="num">Q1 bp</th><th class="num">Median bp</th>
    <th class="num">Q3 bp</th><th class="num">Longest bp</th></tr></thead>
    <tbody>{''.join(rows)}</tbody>
  </table>
</details>
"""
    return chart, table


def _flag_note(depth_flagged: bool, length_flagged: bool) -> str:
    """
    Describe which anomaly flags apply to a contig end.

    Depth and length are scored separately, so an end can carry either, both,
    or neither. Naming the combination is the point: both together is the
    collapsed-array signature, whereas either alone has a more ordinary
    explanation.

    Parameters
    ----------
    depth_flagged : bool
        Whether overhang depth was flagged as anomalous.
    length_flagged : bool
        Whether median overhang length was flagged as anomalous.

    Returns
    -------
    str
        A short description, or an en dash when neither flag applies.
    """
    if depth_flagged and length_flagged:
        return 'anomalous depth and length'
    if depth_flagged:
        return 'anomalous depth'
    if length_flagged:
        return 'anomalous length'
    return '–'


def depth_vs_length_chart(
    stats_dict: Dict[str, ContigStats],
    flagged_left: Sequence[str],
    flagged_right: Sequence[str],
    flagged_left_length: Sequence[str] = (),
    flagged_right_length: Sequence[str] = (),
) -> Tuple[str, str]:
    """
    Plot overhang depth against median overhang length for each contig end.

    The two preceding charts each show one measure against contig order, which
    answers "which end stands out" but not "in what way". Putting the measures
    on the two axes separates the cases that matter: a collapsed repeat or rDNA
    array lands top-right, with many reads each running a long way past the
    end; an organellar contig or a repeat pulling in reads from elsewhere lands
    top-left, deep but short; and a genuine long telomere the assembly stopped
    short of lands bottom-right, which is the case extension exists to serve.

    Both ends of a contig share a ``data-contig`` value, so selecting one point
    highlights its partner here and in the other two charts.

    Parameters
    ----------
    stats_dict : Dict[str, ContigStats]
        Per-contig overhang statistics.
    flagged_left : Sequence[str]
        Contigs flagged for anomalous depth at the left end.
    flagged_right : Sequence[str]
        Contigs flagged for anomalous depth at the right end.
    flagged_left_length : Sequence[str], optional
        Contigs flagged for anomalous overhang length at the left end.
    flagged_right_length : Sequence[str], optional
        Contigs flagged for anomalous overhang length at the right end.

    Returns
    -------
    Tuple[str, str]
        The chart markup and the equivalent table markup. Both are empty when
        no contig end carries an overhang.
    """
    left_flags, right_flags = set(flagged_left), set(flagged_right)
    left_len_flags = set(flagged_left_length)
    right_len_flags = set(flagged_right_length)

    points = []
    for name, stats in stats_dict.items():
        for end, overhangs, depth_flagged, length_flagged in (
            (
                'left',
                stats.left_overhangs,
                name in left_flags,
                name in left_len_flags,
            ),
            (
                'right',
                stats.right_overhangs,
                name in right_flags,
                name in right_len_flags,
            ),
        ):
            if not overhangs:
                continue
            lengths = sorted(float(oh.length) for oh in overhangs)
            points.append(
                {
                    'contig': name,
                    'contig_length': stats.contig_length,
                    'end': end,
                    'depth': len(overhangs),
                    'median': statistics.median(lengths),
                    'longest': lengths[-1],
                    'best_gain': max(oh.net_gain for oh in overhangs),
                    'depth_flagged': depth_flagged,
                    'length_flagged': length_flagged,
                }
            )

    if not points:
        return '', ''

    # Domains are padded past the extremes so that the deepest and longest
    # ends — the ones the chart exists to show — sit inside the plot rather
    # than half on the axis rule.
    max_depth = (max(p['depth'] for p in points) or 1) * 1.08
    max_median = (max(p['median'] for p in points) or 1.0) * 1.08

    plot_w = TARGET_W
    width = plot_w + PAD_L + PAD_R
    height = PLOT_H + PAD_T + PAD_B

    def x_of(median: float) -> float:
        """
        Horizontal position for a median overhang length.

        Parameters
        ----------
        median : float
            Median overhang length in bases.

        Returns
        -------
        float
            X coordinate.
        """
        return PAD_L + (median / max_median) * plot_w

    def y_of(depth: int) -> float:
        """
        Vertical position for an overhang depth.

        Parameters
        ----------
        depth : int
            Number of overhang reads.

        Returns
        -------
        float
            Y coordinate.
        """
        return PAD_T + PLOT_H - (depth / max_depth) * PLOT_H

    svg: List[str] = []

    for frac in (0, 0.25, 0.5, 0.75, 1.0):
        y = PAD_T + PLOT_H - frac * PLOT_H
        svg.append(
            f'<line class="grid" x1="{PAD_L}" y1="{y:.1f}" '
            f'x2="{PAD_L + plot_w}" y2="{y:.1f}"/>'
        )
        svg.append(
            f'<text class="tick" x="{PAD_L - 8}" y="{y + 3.5:.1f}" '
            f'text-anchor="end">{fmt_int(round(max_depth * frac))}</text>'
        )

    # X ticks along the baseline, so the length scale can be read directly.
    for frac in (0.25, 0.5, 0.75, 1.0):
        x = PAD_L + frac * plot_w
        svg.append(
            f'<text class="tick" x="{x:.1f}" y="{PAD_T + PLOT_H + 16}" '
            f'text-anchor="middle">{fmt_int(round(max_median * frac))}</text>'
        )

    for point in points:
        cx, cy = x_of(point['median']), y_of(point['depth'])
        flagged = point['depth_flagged'] or point['length_flagged']

        note = _flag_note(point['depth_flagged'], point['length_flagged'])

        classes = f'pt pt-{point["end"]}' + (' pt-flag' if flagged else '')
        title = (
            f'{escape(point["contig"])} — {point["end"]} end — '
            f'{fmt_int(point["depth"])} reads, '
            f'median {fmt_int(round(point["median"]))} bp'
        )
        svg.append(
            f'<g class="{classes}" tabindex="0" role="listitem"'
            f' data-contig="{escape(point["contig"])}" data-end="{point["end"]}"'
            f' data-depth="{fmt_int(point["depth"])}"'
            f' data-contiglength="{fmt_int(point["contig_length"])}"'
            f' data-median="{fmt_int(round(point["median"]))}"'
            f' data-longest="{fmt_int(round(point["longest"]))}"'
            f' data-gain="{fmt_int(point["best_gain"])}"'
            f' data-flagged="{"yes" if flagged else "no"}"'
            f' data-note="{escape(note)}"'
            f' data-x="{cx:.1f}" data-y="{cy:.1f}">'
            f'<title>{title}</title>'
            f'<circle class="hit" cx="{cx:.1f}" cy="{cy:.1f}" r="12"/>'
            f'<circle class="dot" cx="{cx:.1f}" cy="{cy:.1f}" r="4.5"/>'
            + (
                f'<circle class="ring" cx="{cx:.1f}" cy="{cy:.1f}" r="8.5"/>'
                if flagged
                else ''
            )
            + '</g>'
        )

    for point in points:
        if not (point['depth_flagged'] or point['length_flagged']):
            continue
        name = point['contig']
        half_label = 3.5 * len(name)
        lx = min(
            max(x_of(point['median']), PAD_L + half_label),
            PAD_L + plot_w - half_label,
        )
        svg.append(
            f'<text class="pt-label" x="{lx:.1f}" '
            f'y="{y_of(point["depth"]) - 13:.1f}" '
            f'text-anchor="middle">{escape(name)}</text>'
        )

    svg.append(
        f'<line class="axis" x1="{PAD_L}" y1="{PAD_T + PLOT_H}" '
        f'x2="{PAD_L + plot_w}" y2="{PAD_T + PLOT_H}"/>'
    )
    svg.append(
        f'<line class="axis" x1="{PAD_L}" y1="{PAD_T}" '
        f'x2="{PAD_L}" y2="{PAD_T + PLOT_H}"/>'
    )
    svg.append(
        f'<text class="axis-title" x="{PAD_L}" y="{height - 8}">'
        f'Horizontal: median overhang length (bp). '
        f'Vertical: overhang reads. Ends high on both sit top right.</text>'
    )

    chart = f"""
<div class="chart-scroll">
  <svg viewBox="0 0 {width} {height}" style="min-width:{width}px;width:100%"
       height="{height}" preserveAspectRatio="xMinYMid meet"
       role="img"
       aria-label="Overhang read depth against median overhang length for each contig end">
    {''.join(svg)}
  </svg>
</div>
"""

    rows = ''.join(
        f'<tr{" class=is-flagged" if (p["depth_flagged"] or p["length_flagged"]) else ""}>'
        f'<td>{escape(p["contig"])}</td><td>{p["end"]}</td>'
        f'<td class="num">{fmt_int(p["depth"])}</td>'
        f'<td class="num">{fmt_int(round(p["median"]))}</td>'
        f'<td class="num">{fmt_int(round(p["longest"]))}</td>'
        f'<td class="num">{fmt_int(p["best_gain"])}</td>'
        f'<td>{_flag_note(p["depth_flagged"], p["length_flagged"])}</td>'
        f'</tr>'
        for p in points
    )
    table = f"""
<details class="table-view">
  <summary>Table view</summary>
  <table>
    <thead><tr><th>Contig</th><th>End</th><th class="num">Reads</th>
    <th class="num">Median bp</th><th class="num">Longest bp</th>
    <th class="num">Best gain bp</th><th>Flagged</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
</details>
"""
    return chart, table
