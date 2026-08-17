"""
Stylesheet for the ``teloclip extend`` HTML report.

The palette and the CSS live together because they are one decision: which
colours carry meaning, and what they mean in each mode. Both light and dark
variants are emitted into every report so the file renders correctly wherever
it is opened, with no preference stored and nothing fetched.
"""

from typing import Dict

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


def build_css() -> str:
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

/* Violins. Identity is carried by the outline at full strength; the fill is
   held well back so that two halves meeting at the centre line read as one
   shape rather than as a solid block. */
svg .vio-body {{ stroke-width: 1.5; }}
svg .vio-left .vio-body {{ fill: color-mix(in srgb, var(--left) 22%, transparent);
  stroke: var(--left); }}
svg .vio-right .vio-body {{ fill: color-mix(in srgb, var(--right) 22%, transparent);
  stroke: var(--right); }}
svg .vio-left .vio-dot {{ fill: var(--left); stroke: var(--surface); stroke-width: 2; }}
svg .vio-right .vio-dot {{ fill: var(--right); stroke: var(--surface); stroke-width: 2; }}
/* The interquartile rule and the median tick sit over the density, and have to
   hold their own against the fill beneath them. */
svg .vio-iqr {{ stroke: var(--ink2); stroke-width: 3; stroke-linecap: round; }}
svg .vio-median {{ stroke: var(--ink); stroke-width: 2; stroke-linecap: round; }}
svg .vio-flag .vio-body {{ stroke: var(--critical); }}
svg .vio {{ cursor: pointer; }}
svg .vio:hover .vio-body {{ stroke-width: 2.5; }}

/* Cross-chart selection. Clicking a contig in any chart marks both of its ends
   in all three, so an end that looks unremarkable in one view can be found in
   the others. Everything unselected recedes rather than disappearing: the
   distribution it belongs to is still the context for the selection. */
svg.has-selection .pt:not(.is-selected),
svg.has-selection .vio:not(.is-selected) {{ opacity: .25; }}
svg .pt.is-selected .dot {{ r: 6.5; }}
svg .is-selected .sel-ring {{
  fill: none; stroke: var(--ink); stroke-width: 2;
}}
svg .vio.is-selected .vio-body {{ stroke-width: 2.5; }}

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
.aln-ruler {{ position: relative; height: 2.6em; }}
.aln-ruler .aln-seq {{ position: relative; display: inline-block; height: 100%; }}
.ruler .tick {{
  /* No font-size here: `left` is expressed in `ch`, which must resolve
     against the 12px monospace grid the sequences use. Shrinking the tick
     would shrink its `ch` and drift the scale off the bases.
     The mark occupies only the lower part of the row so the label can sit
     above it with clear air between the numbers and the first read. */
  position: absolute; top: 1.5em; bottom: 0;
  border-left: 1px solid var(--axis);
}}
.ruler .tick i {{
  position: absolute; bottom: 100%; left: 0; padding: 0 0 3px 3px;
  font-style: normal; font-size: 10px; line-height: 1.1;
  color: var(--muted); white-space: nowrap;
}}
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
footer {{
  margin-top: 3rem; padding-top: 1rem; border-top: 1px solid var(--grid);
  color: var(--ink2); font-size: .8rem;
}}
footer dl {{ display: grid; grid-template-columns: max-content 1fr;
  gap: .3rem .9rem; margin: 0; }}
footer dt {{ color: var(--muted); text-transform: uppercase;
  letter-spacing: .04em; font-size: .72rem; padding-top: .1rem; }}
footer dd {{ margin: 0; }}
footer code {{
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-size: .78rem; white-space: pre-wrap; word-break: break-all;
  color: var(--ink);
}}
"""
