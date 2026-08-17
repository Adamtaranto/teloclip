"""
Nature-journal-style figures for the terminal coverage-decay model.

These produce the static images embedded in the documentation
(``docs/guide/coverage-decay.md``). matplotlib is imported lazily inside the
plotting calls so that importing :mod:`teloclip.model` for the analytic maths
never drags a plotting stack in.

Colour anchors on the palette of the ``teloclip extend`` HTML report
(``teloclip.report.css``): blue and orange are its categorical hues, extended
here with a purple and a green (CVD-validated as a set) so each size-selection
cutoff in the family plot gets its own legend-keyed colour.
"""

from contextlib import contextmanager
import logging
from typing import Iterator, Tuple

import numpy as np

from teloclip.model.decay import (
    FixedLength,
    LognormalLengths,
    bootstrap_coverage,
    expected_relative_coverage,
)

logger = logging.getLogger(__name__)

# Categorical anchors shared with the HTML report palette (light mode).
BLUE = '#2a78d6'
ORANGE = '#eb6834'
INK = '#0b0b0b'
MUTED = '#898781'
GRID = '#e1e0d9'

# Categorical colours for the cutoff family: blue and orange from the report
# palette plus a purple and a green stepped to keep every adjacent pair
# distinguishable under colour-vision deficiency (validated worst-pair
# ΔE 15.3 in this order).
CUTOFF_COLORS = ('#2a78d6', '#eb6834', '#8a5cc2', '#2f8f6b')

# Nature-style figure defaults: sans-serif, small type, thin open axes.
# 89 mm is the single-column width; 183 mm the double column.
SINGLE_COL_IN = 89 / 25.4
DOUBLE_COL_IN = 183 / 25.4

NATURE_RC = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 7,
    'axes.titlesize': 8,
    'axes.labelsize': 7,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'axes.linewidth': 0.6,
    'axes.edgecolor': INK,
    'axes.labelcolor': INK,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'xtick.color': INK,
    'ytick.color': INK,
    'lines.linewidth': 1.0,
    'axes.grid': True,
    'grid.color': GRID,
    'grid.linewidth': 0.4,
    'axes.axisbelow': True,
    'legend.frameon': False,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
}


@contextmanager
def nature_style() -> Iterator[None]:
    """
    Apply the Nature-style matplotlib rc parameters within a context.

    Yields
    ------
    None
        Control returns with the rc parameters applied; they are restored
        when the context exits.
    """
    import matplotlib

    with matplotlib.rc_context(NATURE_RC):
        yield


def plot_decay_families(out_path: str) -> None:
    """
    Plot analytic decay curves across cutoffs and length distributions.

    Two panels: (a) a lognormal read-length profile under increasingly harsh
    size-selection cutoffs, and (b) different fragment-length distributions
    under one fixed cutoff, showing how distribution shape moves the decay.

    Parameters
    ----------
    out_path : str
        Destination PNG path.
    """
    import matplotlib.pyplot as plt

    profile = LognormalLengths(mean=8000.0, sd=4000.0)
    x = np.linspace(0.0, 30_000.0, 600)
    xk = x / 1000.0

    with nature_style():
        fig, (ax_a, ax_b) = plt.subplots(
            1, 2, figsize=(DOUBLE_COL_IN, 2.4), sharey=True
        )

        # Panel a: one distribution, one distinctly coloured curve per
        # cutoff, identified through the legend.
        cutoffs = (0.0, 5000.0, 10_000.0, 20_000.0)
        for cutoff, color in zip(cutoffs, CUTOFF_COLORS):
            c = expected_relative_coverage(x, profile, min_length=cutoff)
            label = 'no cutoff' if cutoff == 0 else f'cutoff ≥{cutoff / 1000:g} kb'
            ax_a.plot(xk, c, color=color, label=label)
        ax_a.legend(loc='lower right')
        ax_a.set_xlabel('Distance from contig end (kb)')
        ax_a.set_ylabel('Relative coverage $c(x)$')
        ax_a.set_title(
            'a  Size-selection cutoff (lognormal 8 ± 4 kb)',
            loc='left',
            fontweight='bold',
        )

        # Panel b: one cutoff, distribution identity in the two categorical
        # hues plus line style, each directly labelled.
        cutoff = 5000.0
        series = (
            (LognormalLengths(mean=8000.0, sd=2000.0), BLUE, '-', 'lognormal, sd 2 kb'),
            (
                LognormalLengths(mean=8000.0, sd=6000.0),
                BLUE,
                '--',
                'lognormal, sd 6 kb',
            ),
            (FixedLength(8000.0), ORANGE, '-', 'fixed 8 kb'),
        )
        for dist, color, style, label in series:
            c = expected_relative_coverage(x, dist, min_length=cutoff)
            ax_b.plot(xk, c, color=color, linestyle=style, label=label)
        ax_b.set_xlabel('Distance from contig end (kb)')
        ax_b.set_title(
            'b  Fragment lengths (cutoff 5 kb)', loc='left', fontweight='bold'
        )
        ax_b.legend(loc='lower right')

        for ax in (ax_a, ax_b):
            ax.set_xlim(0, 30)
            ax.set_ylim(0, 1.05)

        fig.savefig(out_path)
        plt.close(fig)
    logger.info('Wrote decay-family figure to %s', out_path)


def plot_decay_with_ci(
    out_path: str,
    depth: float = 30.0,
    n_reps: int = 200,
    seed: int = 42,
) -> None:
    """
    Plot simulated terminal coverage with a bootstrap band and analytic curve.

    Parameters
    ----------
    out_path : str
        Destination PNG path.
    depth : float, optional
        Interior sequencing depth for the simulation.
    n_reps : int, optional
        Number of bootstrap replicates behind the confidence band.
    seed : int, optional
        Seed for the replicate generators.
    """
    import matplotlib.pyplot as plt

    profile = LognormalLengths(mean=8000.0, sd=4000.0)
    cutoff = 2500.0
    window = 30_000
    result = bootstrap_coverage(
        chrom_length=1_000_000,
        depth=depth,
        dist=profile,
        min_length=cutoff,
        window=window,
        n_reps=n_reps,
        seed=seed,
    )
    xk = result.x / 1000.0

    with nature_style():
        fig, ax = plt.subplots(figsize=(SINGLE_COL_IN, 2.3))

        ax.fill_between(
            xk,
            result.lo,
            result.hi,
            color=BLUE,
            alpha=0.18,
            linewidth=0,
            label=f'95% band ({n_reps} simulations)',
        )
        ax.plot(xk, result.mean, color=BLUE, linewidth=0.8, label='simulated mean')
        ax.plot(
            xk,
            result.analytic,
            color=ORANGE,
            linestyle='--',
            linewidth=1.0,
            label='analytic $D\\,c(x)$',
        )
        ax.axhline(depth, color=MUTED, linewidth=0.5, linestyle=':')
        ax.annotate(
            f'interior depth {depth:g}×',
            (0.98, depth),
            xycoords=('axes fraction', 'data'),
            xytext=(0, 3),
            textcoords='offset points',
            ha='right',
            fontsize=6,
            color=MUTED,
        )

        ax.set_xlabel('Distance from contig end (kb)')
        ax.set_ylabel('Coverage (reads)')
        ax.set_xlim(0, window / 1000.0)
        ax.set_ylim(0, None)
        ax.set_title(
            f'Lognormal fragments {profile.mean / 1000:g} ± {profile.sd / 1000:g} kb'
            f' · cutoff ≥{cutoff / 1000:g} kb · depth {depth:g}×',
            loc='left',
            fontweight='bold',
        )
        ax.legend(loc='lower right')

        fig.savefig(out_path)
        plt.close(fig)
    logger.info('Wrote bootstrap-CI figure to %s', out_path)


def generate_all(out_dir: str) -> Tuple[str, str]:
    """
    Generate every documentation figure into a directory.

    Parameters
    ----------
    out_dir : str
        Directory to write the PNGs into (must exist).

    Returns
    -------
    Tuple[str, str]
        Paths of the decay-family figure and the bootstrap-CI figure.
    """
    import os

    families = os.path.join(out_dir, 'coverage-decay-analytic.png')
    ci = os.path.join(out_dir, 'coverage-decay-bootstrap.png')
    plot_decay_families(families)
    plot_decay_with_ci(ci)
    return families, ci
