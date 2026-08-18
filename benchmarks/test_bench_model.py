"""
Benchmarks for the terminal coverage-decay model.

The analytic curve is evaluated per point in the docs figures and tests, and
the Monte-Carlo simulator dominates the bootstrap figure's cost, so both are
tracked. Sizes mirror the committed documentation figure (megabase molecule,
tens of thousands of fragments).
"""

import pytest

np = pytest.importorskip('numpy')

from teloclip.model.decay import (  # noqa: E402
    LognormalLengths,
    bootstrap_coverage,
    expected_relative_coverage,
    simulate_coverage,
)

PROFILE = LognormalLengths(mean=8000.0, sd=4000.0)


@pytest.mark.benchmark
def test_analytic_curve_dense_grid():
    """Benchmark the closed-form curve on a 10k-point grid."""
    x = np.linspace(0.0, 100_000.0, 10_000)
    for _ in range(10):
        expected_relative_coverage(x, PROFILE, min_length=5000.0)


@pytest.mark.benchmark
def test_simulate_megabase_molecule():
    """Benchmark one simulation with ~1e5 fragments (5 Mb at 160x depth)."""
    rng = np.random.default_rng(0)
    simulate_coverage(
        chrom_length=5_000_000,
        depth=160.0,
        dist=PROFILE,
        min_length=2500.0,
        window=50_000,
        rng=rng,
    )


@pytest.mark.benchmark
def test_bootstrap_small():
    """Benchmark a small bootstrap: ten replicates with their percentile band."""
    bootstrap_coverage(
        chrom_length=1_000_000,
        depth=30.0,
        dist=PROFILE,
        min_length=2500.0,
        window=30_000,
        n_reps=10,
        seed=7,
    )
