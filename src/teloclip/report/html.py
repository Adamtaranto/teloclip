"""
Self-contained HTML report for ``teloclip extend``.

The Markdown report answers "what changed". This one answers "should I believe
it", by showing the evidence: every overhang read laid out against the contig
terminus it supports, and the distribution of overhang depth and length across
the assembly so anomalies are visible rather than merely listed.

The output is a single file with no external references — no CDN scripts, no
web fonts, no separate stylesheet — so it survives being emailed, copied onto a
cluster, or opened years later.

This module assembles the document. The pieces come from
:mod:`teloclip.report.panels` (per-read alignment panels),
:mod:`teloclip.report.charts` (inline SVG charts) and
:mod:`teloclip.report.css` (the stylesheet).
"""

from html import escape
from typing import Dict, List, Sequence

from ..core.analysis import ContigStats
from .charts import coverage_chart
from .css import build_css
from .panels import motif_pattern
from .text import fmt_delta_html, fmt_int


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
    generated: str = '',
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
        Command line that produced the report, verbatim.
    generated : str, optional
        Timestamp for the footer.

    Returns
    -------
    str
        A complete, self-contained HTML document.
    """
    pattern = motif_pattern(motifs)

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
        ('Contigs', fmt_int(total_contigs)),
        ('Contigs extended', fmt_int(len(extensions_applied))),
        ('Ends extended', f'{fmt_int(ends_extended)} of {fmt_int(total_contigs * 2)}'),
        ('Net bases gained', fmt_int(net_gained)),
        ('Bases trimmed', fmt_int(trimmed)),
        ('Flagged contigs', fmt_int(len(flagged_all))),
        ('Warnings', fmt_int(len(warnings))),
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
            f'<td class="num">{fmt_int(info.get("original_length", 0))}</td>'
            f'<td class="num">{fmt_int(info.get("final_length", 0))}</td>'
            f'<td class="num">{fmt_delta_html(left)}</td>'
            f'<td class="num">{fmt_delta_html(right)}</td>'
            f'<td class="num">{fmt_delta_html(left + right)}</td>'
            f'<td class="num">{fmt_int(info.get("left_trim_length", 0))}</td>'
            f'<td class="num">{fmt_int(info.get("right_trim_length", 0))}</td></tr>'
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
    chart, chart_table = coverage_chart(
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

    # Provenance: what produced this file, and exactly how. Labelled rather
    # than run together, so the version cannot be mistaken for part of the
    # command line.
    footer_rows = [('Version', f'teloclip {version}' if version else 'teloclip')]
    footer_rows.append(('Generated', generated))
    if command:
        footer_rows.append(('Command', command))

    footer = (
        '<dl>'
        + ''.join(
            f'<dt>{escape(k)}</dt><dd>'
            + (f'<code>{escape(v)}</code>' if k == 'Command' else escape(v))
            + '</dd>'
            for k, v in footer_rows
        )
        + '</dl>'
    )

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Teloclip extend report</title>
<style>{build_css()}</style>
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
  {'Highlighted bases are motif matches.' if pattern else ''}</p>
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
