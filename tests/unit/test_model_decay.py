"""
Tests for the terminal coverage-decay model.

The analytic curve has strong structural guarantees (exact ramp for fixed
lengths, monotonicity, convergence to interior coverage) and two independent
implementations (closed form and quadrature), so most tests cross-check those
against each other and against the Monte-Carlo simulator.
"""

import pytest

np = pytest.importorskip('numpy')

from teloclip.model import (  # noqa: E402
    BootstrapResult,
    FixedLength,
    LognormalLengths,
    bootstrap_coverage,
    expected_relative_coverage,
    expected_relative_coverage_quadrature,
    simulate_coverage,
)

# A realistic long-read profile: 8 kb mean, 4 kb sd, in bases.
REALISTIC = LognormalLengths(mean=8000.0, sd=4000.0)


class TestDistributions:
    """Construction and moments of the length distributions."""

    def test_fixed_rejects_nonpositive_length(self):
        """A degenerate distribution needs a positive length."""
        with pytest.raises(ValueError, match='positive'):
            FixedLength(0)

    @pytest.mark.parametrize('mean,sd', [(0, 100), (5000, 0), (-1, 1)])
    def test_lognormal_rejects_nonpositive_params(self, mean, sd):
        """Real-space mean and sd must both be positive."""
        with pytest.raises(ValueError, match='positive'):
            LognormalLengths(mean=mean, sd=sd)

    def test_lognormal_real_space_parameterisation(self):
        """log_mu/log_sigma must reproduce the requested real-space moments."""
        dist = REALISTIC
        mean = np.exp(dist.log_mu + dist.log_sigma**2 / 2)
        var = (np.exp(dist.log_sigma**2) - 1) * mean**2
        assert mean == pytest.approx(8000.0)
        assert np.sqrt(var) == pytest.approx(4000.0)

    def test_lognormal_samples_match_moments(self):
        """Large seeded samples should land near the requested moments."""
        rng = np.random.default_rng(7)
        lens = REALISTIC.sample(rng, 200_000)
        assert lens.mean() == pytest.approx(8000.0, rel=0.02)
        assert lens.std() == pytest.approx(4000.0, rel=0.05)

    def test_survival_is_monotone_and_bounded(self):
        """S(u) must fall from 1 towards 0 and stay in [0, 1]."""
        u = np.linspace(0, 100_000, 500)
        s = REALISTIC.survival(u)
        assert np.all(np.diff(s) <= 0)
        assert s[0] == pytest.approx(1.0)
        assert 0.0 <= s[-1] < 1e-6  # ~12 mean lengths into the tail

    def test_truncated_mean_untruncated_matches_mean(self):
        """A zero cutoff must leave the mean untouched."""
        assert REALISTIC.truncated_mean(0.0) == pytest.approx(8000.0)
        assert FixedLength(5000).truncated_mean(0.0) == 5000

    def test_truncated_mean_increases_with_cutoff(self):
        """Removing short fragments can only raise the mean length."""
        cutoffs = [0, 1000, 4000, 8000, 15000]
        means = [REALISTIC.truncated_mean(c) for c in cutoffs]
        assert np.all(np.diff(means) > 0)
        # And always at least the cutoff itself.
        for cutoff, mean in zip(cutoffs, means):
            assert mean >= max(cutoff, 8000.0)

    def test_truncated_mean_extreme_cutoff_raises(self):
        """A cutoff beyond the whole distribution has no surviving fragments."""
        with pytest.raises(ValueError, match='survive'):
            REALISTIC.truncated_mean(10_000_000.0)
        with pytest.raises(ValueError, match='survive'):
            FixedLength(5000).truncated_mean(5001)


class TestAnalyticCurve:
    """Closed-form expected relative coverage."""

    def test_fixed_length_is_exact_ramp(self):
        """c(x) = min(x + 1, L) / L exactly for a fixed length."""
        length = 5000
        x = np.arange(0, 12_000, 7)
        c = expected_relative_coverage(x, FixedLength(length))
        np.testing.assert_allclose(c, np.minimum(x + 1, length) / length, rtol=0)

    def test_fixed_length_ignores_permissive_cutoff(self):
        """Any cutoff below L keeps every fragment: the ramp is unchanged."""
        x = np.arange(0, 8000, 11)
        base = expected_relative_coverage(x, FixedLength(5000))
        cut = expected_relative_coverage(x, FixedLength(5000), min_length=3000)
        np.testing.assert_array_equal(base, cut)

    def test_fixed_length_impossible_cutoff_raises(self):
        """A cutoff above L discards everything."""
        with pytest.raises(ValueError, match='survive'):
            expected_relative_coverage(
                np.arange(10.0), FixedLength(5000), min_length=6000
            )

    @pytest.mark.parametrize('min_length', [0.0, 1000.0, 8000.0])
    def test_monotone_and_bounded(self, min_length):
        """c(x) rises from ~0 towards 1 and never leaves [0, 1]."""
        x = np.linspace(0, 100_000, 2000)
        c = expected_relative_coverage(x, REALISTIC, min_length=min_length)
        assert np.all(np.diff(c) >= 0)
        assert np.all(c >= 0.0)
        assert np.all(c <= 1.0 + 1e-12)

    def test_converges_to_interior_coverage(self):
        """Far from the terminus the relative coverage must reach 1."""
        c = expected_relative_coverage(np.array([500_000.0]), REALISTIC)
        assert c[0] == pytest.approx(1.0, abs=1e-9)

    def test_coverage_at_terminus(self):
        """The outermost base is covered only by fragments starting on it."""
        c = expected_relative_coverage(np.array([0.0]), REALISTIC)
        assert c[0] == pytest.approx(1.0 / 8000.0, rel=1e-6)

    def test_zero_cutoff_matches_untruncated_branch(self):
        """min_length=0 must reduce to the plain partial-expectation formula."""
        x = np.linspace(0, 40_000, 500)
        a = expected_relative_coverage(x, REALISTIC, min_length=0.0)
        b = expected_relative_coverage(x, REALISTIC, min_length=1e-9)
        np.testing.assert_allclose(a, b, rtol=1e-9)

    def test_size_selection_deepens_the_decay(self):
        """A harsher cutoff must depress terminal coverage further."""
        x = np.linspace(0, 6000, 300)
        mild = expected_relative_coverage(x, REALISTIC, min_length=1000)
        harsh = expected_relative_coverage(x, REALISTIC, min_length=10_000)
        assert np.all(harsh <= mild + 1e-12)
        # Strictly worse close to the terminus.
        assert harsh[10] < mild[10]

    @pytest.mark.parametrize(
        'dist,min_length',
        [
            (REALISTIC, 0.0),
            (REALISTIC, 2000.0),
            (REALISTIC, 12_000.0),
            (LognormalLengths(mean=15_000, sd=12_000), 5000.0),
            (LognormalLengths(mean=500, sd=100), 0.0),
        ],
    )
    def test_closed_form_matches_quadrature(self, dist, min_length):
        """The two independent implementations must agree closely."""
        x = np.linspace(0, 12 * dist.mean, 400)
        closed = expected_relative_coverage(x, dist, min_length=min_length)
        quad = expected_relative_coverage_quadrature(
            x, dist.survival, min_length=min_length, n_nodes=16_384
        )
        np.testing.assert_allclose(quad, closed, rtol=2e-5, atol=1e-7)

    def test_quadrature_handles_fixed_length_survival(self):
        """The generic path also reproduces the fixed-length ramp."""
        length = 5000.0
        x = np.linspace(0, 9000, 300)
        quad = expected_relative_coverage_quadrature(
            x,
            FixedLength(length).survival,
            n_nodes=32_768,
        )
        np.testing.assert_allclose(quad, np.minimum(x + 1, length) / length, atol=5e-4)

    def test_quadrature_extreme_cutoff_raises(self):
        """The generic path applies the same no-survivors guard."""
        with pytest.raises(ValueError, match='survive'):
            expected_relative_coverage_quadrature(
                np.arange(10.0), REALISTIC.survival, min_length=10_000_000.0
            )


def _mean_coverage(seed: int, n_reps: int, **kwargs) -> 'np.ndarray':
    """
    Average several independent simulation replicates.

    Coverage is correlated over the fragment length (~kb), so a single
    replicate has large regional fluctuations even at high depth; averaging
    replicates shrinks them by sqrt(n_reps) and lets tests use tight bounds.

    Parameters
    ----------
    seed : int
        Seed for the spawned replicate generators.
    n_reps : int
        Number of independent replicates to average.
    **kwargs
        Forwarded to :func:`simulate_coverage`.

    Returns
    -------
    np.ndarray
        Per-position mean coverage across replicates.
    """
    root = np.random.default_rng(seed)
    reps = [simulate_coverage(rng=child, **kwargs) for child in root.spawn(n_reps)]
    return np.mean(reps, axis=0)


class TestSimulator:
    """Monte-Carlo simulation of terminal coverage."""

    def test_seed_determinism(self):
        """Identical seeds must give identical coverage arrays."""
        kwargs = {
            'chrom_length': 1_000_000,
            'depth': 20.0,
            'dist': REALISTIC,
            'min_length': 1000.0,
            'window': 20_000,
        }
        a = simulate_coverage(rng=np.random.default_rng(11), **kwargs)
        b = simulate_coverage(rng=np.random.default_rng(11), **kwargs)
        np.testing.assert_array_equal(a, b)

    def test_window_larger_than_molecule_raises(self):
        """The tallied window has to fit inside the molecule."""
        with pytest.raises(ValueError, match='exceeds'):
            simulate_coverage(
                chrom_length=10_000,
                depth=5.0,
                dist=REALISTIC,
                min_length=0.0,
                window=20_000,
                rng=np.random.default_rng(0),
            )

    def test_interior_depth_matches_target(self):
        """Deep in the window the mean coverage should approach depth D."""
        cov = _mean_coverage(
            seed=3,
            n_reps=30,
            chrom_length=2_000_000,
            depth=30.0,
            dist=REALISTIC,
            min_length=1000.0,
            window=100_000,
        )
        interior = cov[60_000:]  # far past the decay region (~7 mean lengths)
        # Per-replicate regional means have sd ~2.5, so 30 replicates put the
        # standard error near 0.46; a 7% bound sits beyond 4 sigma.
        assert interior.mean() == pytest.approx(30.0, rel=0.07)

    def test_matches_analytic_curve(self):
        """A deep simulation should track D * c(x) closely."""
        depth = 200.0
        window = 40_000
        cov = _mean_coverage(
            seed=5,
            n_reps=20,
            chrom_length=2_000_000,
            depth=depth,
            dist=REALISTIC,
            min_length=2000.0,
            window=window,
        )
        x = np.arange(window, dtype=float)
        expected = depth * expected_relative_coverage(x, REALISTIC, min_length=2000.0)
        # Compare in 4 kb bins to average out fragment-scale noise.
        bins = 10
        sim_binned = cov.reshape(bins, -1).mean(axis=1)
        exp_binned = expected.reshape(bins, -1).mean(axis=1)
        np.testing.assert_allclose(sim_binned, exp_binned, rtol=0.06, atol=1.0)

    def test_fixed_length_terminal_ramp(self):
        """With fixed lengths the simulated decay is the linear ramp."""
        length = 5000
        depth = 300.0
        window = 8000
        cov = _mean_coverage(
            seed=9,
            n_reps=25,
            chrom_length=1_000_000,
            depth=depth,
            dist=FixedLength(length),
            min_length=0.0,
            window=window,
        )
        x = np.arange(window, dtype=float)
        expected = depth * np.minimum(x + 1, length) / length
        bins = 16  # 500 bp bins
        np.testing.assert_allclose(
            cov.reshape(bins, -1).mean(axis=1),
            expected.reshape(bins, -1).mean(axis=1),
            rtol=0.06,
            atol=2.0,
        )

    def test_coverage_depleted_at_terminus(self):
        """The first base must be far below interior depth."""
        cov = simulate_coverage(
            chrom_length=1_000_000,
            depth=50.0,
            dist=REALISTIC,
            min_length=5000.0,
            window=30_000,
            rng=np.random.default_rng(13),
        )
        assert cov[0] < 5
        assert cov[:100].mean() < 0.2 * cov[-5000:].mean()


@pytest.fixture(scope='module')
def result() -> 'BootstrapResult':
    """One moderately sized bootstrap shared across the band assertions."""
    return bootstrap_coverage(
        chrom_length=500_000,
        depth=25.0,
        dist=REALISTIC,
        min_length=1000.0,
        window=15_000,
        n_reps=25,
        seed=42,
    )


class TestBootstrap:
    """Bootstrap confidence bands around the simulated decay."""

    def test_shapes_are_consistent(self, result):
        """Every band component covers the same positions."""
        n = result.x.size
        assert (
            result.mean.size
            == result.lo.size
            == result.hi.size
            == result.analytic.size
            == n
        )

    def test_band_orders_correctly(self, result):
        """Per position: lo <= mean <= hi."""
        assert np.all(result.lo <= result.mean + 1e-12)
        assert np.all(result.mean <= result.hi + 1e-12)

    def test_band_contains_analytic_curve(self, result):
        """The analytic expectation should sit inside the band almost always."""
        inside = (result.analytic >= result.lo) & (result.analytic <= result.hi)
        assert inside.mean() > 0.9

    def test_seed_determinism(self, result):
        """The same seed must reproduce the same band."""
        again = bootstrap_coverage(
            chrom_length=500_000,
            depth=25.0,
            dist=REALISTIC,
            min_length=1000.0,
            window=15_000,
            n_reps=25,
            seed=42,
        )
        np.testing.assert_array_equal(result.mean, again.mean)
        np.testing.assert_array_equal(result.lo, again.lo)
        np.testing.assert_array_equal(result.hi, again.hi)
