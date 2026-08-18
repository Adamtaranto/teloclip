"""
Expected coverage decay towards the end of a linear DNA molecule.

Random fragmentation places fragment start points uniformly along a molecule,
so a base at distance ``x`` from a terminus can only be covered by fragments
that start within the last ``x`` bases before it — at most ``min(x + 1, l)``
placements for a fragment of length ``l``, against ``l`` placements for an
interior base. Size selection then discards fragments shorter than a cutoff,
which removes exactly the fragments most likely to cover terminal bases.
Averaging over the (truncated) fragment-length distribution gives the expected
coverage at ``x`` relative to the interior:

    c(x) = E[min(x + 1, l) | l >= Lmin] / E[l | l >= Lmin]

Everything here works with the survival function ``S(u) = P(l > u)`` through
the identity ``E[min(a, l)] = integral of S(u) for u in [0, a]``, which keeps
the maths closed-form for the fixed-length and lognormal cases and reduces the
general case to one numerical integral.

The module deliberately assumes a semi-infinite molecule (one terminus, with
the interior far longer than any fragment): both-ends effects on a real
chromosome are just this same decay mirrored at the other terminus.
"""

from dataclasses import dataclass
import logging
import math
from typing import Callable, Union

import numpy as np

logger = logging.getLogger(__name__)

# math.erf vectorised over arrays. numpy has no erf of its own and scipy is
# not a dependency; the per-element Python call is irrelevant next to the
# simulation costs in this module.
_erf = np.vectorize(math.erf, otypes=[float])

# Refuse cutoffs that discard essentially the whole length distribution:
# below this survival probability the truncated moments are numerically
# meaningless (and the experiment they describe would yield no library).
_MIN_SURVIVAL = 1e-12


def _norm_cdf(z: np.ndarray) -> np.ndarray:
    """
    Evaluate the standard normal cumulative distribution function.

    Parameters
    ----------
    z : np.ndarray
        Points at which to evaluate the CDF.

    Returns
    -------
    np.ndarray
        ``P(Z <= z)`` for a standard normal ``Z``.
    """
    return 0.5 * (1.0 + _erf(np.asarray(z, dtype=float) / math.sqrt(2.0)))


@dataclass(frozen=True)
class FixedLength:
    """
    Degenerate fragment-length distribution: every fragment has one length.

    Attributes
    ----------
    length : float
        Fragment length in bases. Must be positive.
    """

    length: float

    def __post_init__(self) -> None:
        """Validate the fragment length."""
        if self.length <= 0:
            raise ValueError(f'Fragment length must be positive, got {self.length}')

    @property
    def mean(self) -> float:
        """
        Mean fragment length.

        Returns
        -------
        float
            The fixed fragment length.
        """
        return self.length

    def survival(self, u: np.ndarray) -> np.ndarray:
        """
        Survival function ``P(l > u)``.

        Parameters
        ----------
        u : np.ndarray
            Lengths at which to evaluate the survival function.

        Returns
        -------
        np.ndarray
            1.0 where ``u < length``, else 0.0.
        """
        return (np.asarray(u, dtype=float) < self.length).astype(float)

    def truncated_mean(self, min_length: float) -> float:
        """
        Mean fragment length after size selection at ``min_length``.

        Parameters
        ----------
        min_length : float
            Size-selection cutoff; fragments shorter than this are discarded.

        Returns
        -------
        float
            The fixed length (truncation keeps every fragment or none).

        Raises
        ------
        ValueError
            If the cutoff exceeds the fragment length, so that no fragment
            survives size selection.
        """
        if min_length > self.length:
            raise ValueError(
                f'Size-selection cutoff {min_length} exceeds the fixed fragment '
                f'length {self.length}: no fragments survive'
            )
        return self.length

    def sample(self, rng: np.random.Generator, n: int) -> np.ndarray:
        """
        Draw fragment lengths.

        Parameters
        ----------
        rng : np.random.Generator
            Source of randomness (unused for the degenerate distribution,
            accepted for interface symmetry).
        n : int
            Number of fragments to draw.

        Returns
        -------
        np.ndarray
            Array of ``n`` copies of the fixed length.
        """
        return np.full(n, self.length, dtype=float)


@dataclass(frozen=True)
class LognormalLengths:
    """
    Lognormal fragment-length distribution, parameterised in real space.

    The lognormal is the standard empirical model for long-read and sheared
    fragment lengths, and its partial expectation ``E[min(a, l)]`` has a
    closed form in the normal CDF, so the coverage-decay curve needs no
    numerical integration.

    Attributes
    ----------
    mean : float
        Mean fragment length in bases (real space, not log space).
    sd : float
        Standard deviation of fragment length in bases (real space).
    """

    mean: float
    sd: float

    def __post_init__(self) -> None:
        """Validate the real-space parameters."""
        if self.mean <= 0:
            raise ValueError(f'Mean fragment length must be positive, got {self.mean}')
        if self.sd <= 0:
            raise ValueError(f'Fragment length sd must be positive, got {self.sd}')

    @property
    def log_sigma(self) -> float:
        """
        Log-space shape parameter sigma.

        Returns
        -------
        float
            ``sqrt(ln(1 + sd^2 / mean^2))``.
        """
        return math.sqrt(math.log1p((self.sd / self.mean) ** 2))

    @property
    def log_mu(self) -> float:
        """
        Log-space location parameter mu.

        Returns
        -------
        float
            ``ln(mean) - sigma^2 / 2``.
        """
        return math.log(self.mean) - 0.5 * self.log_sigma**2

    def survival(self, u: np.ndarray) -> np.ndarray:
        """
        Survival function ``P(l > u)``.

        Parameters
        ----------
        u : np.ndarray
            Lengths at which to evaluate the survival function.

        Returns
        -------
        np.ndarray
            ``1 - Phi((ln u - mu) / sigma)``, with the value 1.0 for
            ``u <= 0`` where the log is undefined.
        """
        u = np.asarray(u, dtype=float)
        out = np.ones_like(u)
        pos = u > 0
        z = (np.log(u[pos]) - self.log_mu) / self.log_sigma
        out[pos] = 1.0 - _norm_cdf(z)
        return out

    def partial_expectation(self, a: np.ndarray) -> np.ndarray:
        """
        Partial expectation ``E[min(a, l)]``.

        Uses the closed form for the lognormal partial mean:
        ``E[min(a, l)] = m * Phi((ln a - mu - sigma^2) / sigma)
        + a * (1 - Phi((ln a - mu) / sigma))`` where ``m = exp(mu + sigma^2/2)``
        is the distribution mean.

        Parameters
        ----------
        a : np.ndarray
            Capping values; entries at or below zero yield 0.

        Returns
        -------
        np.ndarray
            ``E[min(a, l)]`` for each entry of ``a``.
        """
        a = np.asarray(a, dtype=float)
        out = np.zeros_like(a)
        pos = a > 0
        ap = a[pos]
        z = (np.log(ap) - self.log_mu) / self.log_sigma
        out[pos] = self.mean * _norm_cdf(z - self.log_sigma) + ap * (1.0 - _norm_cdf(z))
        return out

    def truncated_mean(self, min_length: float) -> float:
        """
        Mean fragment length after size selection at ``min_length``.

        Computed as ``Lmin + (E[l] - E[min(Lmin, l)]) / S(Lmin)``, which is
        the integral of the truncated survival function.

        Parameters
        ----------
        min_length : float
            Size-selection cutoff; fragments shorter than this are discarded.

        Returns
        -------
        float
            ``E[l | l >= min_length]``.

        Raises
        ------
        ValueError
            If the cutoff sits so far into the upper tail that essentially
            no fragment survives size selection.
        """
        if min_length <= 0:
            return self.mean
        surv = float(self.survival(np.array([min_length]))[0])
        if surv < _MIN_SURVIVAL:
            raise ValueError(
                f'Size-selection cutoff {min_length} leaves a survival '
                f'probability of {surv:.3g}: effectively no fragments survive'
            )
        capped = float(self.partial_expectation(np.array([min_length]))[0])
        return min_length + (self.mean - capped) / surv

    def sample(self, rng: np.random.Generator, n: int) -> np.ndarray:
        """
        Draw fragment lengths.

        Parameters
        ----------
        rng : np.random.Generator
            Source of randomness.
        n : int
            Number of fragments to draw.

        Returns
        -------
        np.ndarray
            ``n`` lognormal fragment lengths in bases.
        """
        return rng.lognormal(self.log_mu, self.log_sigma, size=n)


LengthDistribution = Union[FixedLength, LognormalLengths]


def expected_relative_coverage(
    x: np.ndarray,
    dist: LengthDistribution,
    min_length: float = 0.0,
) -> np.ndarray:
    """
    Compute expected coverage at distance ``x`` from a terminus.

    Evaluates ``c(x) = E[min(x + 1, l) | l >= Lmin] / E[l | l >= Lmin]``
    in closed form. Multiply by the interior sequencing depth ``D`` to get
    absolute expected coverage.

    Parameters
    ----------
    x : np.ndarray
        Distances from the terminus in bases (0 = outermost base).
    dist : FixedLength or LognormalLengths
        Fragment-length distribution before size selection.
    min_length : float, optional
        Size-selection cutoff: fragments shorter than this are discarded.
        Zero (default) means no size selection.

    Returns
    -------
    np.ndarray
        Relative coverage ``c(x)`` in [0, 1] for each entry of ``x``.

    Raises
    ------
    ValueError
        If ``min_length`` discards (essentially) every fragment.
    """
    x = np.asarray(x, dtype=float)
    a = x + 1.0  # a fragment covering position x must start in the last a bases

    if isinstance(dist, FixedLength):
        # Degenerate case: truncation keeps every fragment or none, so the
        # curve is the plain linear ramp min(a, L) / L.
        dist.truncated_mean(min_length)  # raises if the cutoff kills everything
        return np.minimum(a, dist.length) / dist.length

    denom = dist.truncated_mean(min_length)
    if min_length <= 0:
        return dist.partial_expectation(a) / denom

    # Truncated survival integral: below the cutoff every surviving fragment
    # is longer than u, so the integrand is 1; above it the survival function
    # is renormalised by S(Lmin).
    surv_min = float(dist.survival(np.array([min_length]))[0])
    cap_min = float(dist.partial_expectation(np.array([min_length]))[0])
    numer = np.where(
        a <= min_length,
        a,
        min_length + (dist.partial_expectation(a) - cap_min) / surv_min,
    )
    return numer / denom


def expected_relative_coverage_quadrature(
    x: np.ndarray,
    survival: Callable[[np.ndarray], np.ndarray],
    min_length: float = 0.0,
    upper: Union[float, None] = None,
    n_nodes: int = 8192,
) -> np.ndarray:
    """
    Relative coverage from an arbitrary survival function, by quadrature.

    Generic counterpart to :func:`expected_relative_coverage` for length
    distributions without a closed-form partial expectation (empirical read
    length distributions, gamma, mixtures). Integrates the truncated survival
    function with the trapezoid rule on a grid that includes every query
    point, so no interpolation error is added on top of the quadrature error.

    Parameters
    ----------
    x : np.ndarray
        Distances from the terminus in bases (0 = outermost base).
    survival : Callable[[np.ndarray], np.ndarray]
        Vectorised survival function ``S(u) = P(l > u)`` of the untruncated
        fragment-length distribution.
    min_length : float, optional
        Size-selection cutoff: fragments shorter than this are discarded.
    upper : float or None, optional
        Upper integration limit. Defaults to doubling from the largest query
        point until ``S(upper) < 1e-9``.
    n_nodes : int, optional
        Number of evenly spaced base nodes in the integration grid.

    Returns
    -------
    np.ndarray
        Relative coverage ``c(x)`` for each entry of ``x``.

    Raises
    ------
    ValueError
        If ``min_length`` discards (essentially) every fragment.
    """
    x = np.asarray(x, dtype=float)
    a = x + 1.0

    if upper is None:
        # Walk the tail out until the distribution is exhausted; the trapezoid
        # tail beyond this point contributes < 1e-9 * upper to the integrals.
        upper = max(float(np.max(a)), min_length, 1.0)
        for _ in range(64):
            if float(survival(np.array([upper]))[0]) < 1e-9:
                break
            upper *= 2.0
        else:  # pragma: no cover - only reachable for a defective survival fn
            raise ValueError('Survival function does not decay: cannot bound tail')

    surv_min = float(survival(np.array([min_length]))[0]) if min_length > 0 else 1.0
    if surv_min < _MIN_SURVIVAL:
        raise ValueError(
            f'Size-selection cutoff {min_length} leaves a survival probability '
            f'of {surv_min:.3g}: effectively no fragments survive'
        )

    # One shared grid holding the regular nodes plus every query point and the
    # cutoff itself, so each requested integral ends exactly on a node.
    grid = np.unique(
        np.concatenate(
            [np.linspace(0.0, float(upper), n_nodes), a[a <= upper], [min_length]]
        )
    )
    s_grid = np.asarray(survival(grid), dtype=float)
    # Truncated survival: certain below the cutoff, renormalised above it.
    s_trunc = np.where(grid < min_length, 1.0, s_grid / surv_min)

    # Cumulative trapezoid of the truncated survival function along the grid.
    steps = np.diff(grid)
    increments = steps * 0.5 * (s_trunc[1:] + s_trunc[:-1])
    cumulative = np.concatenate([[0.0], np.cumsum(increments)])

    denom = float(cumulative[-1])
    numer = np.interp(np.minimum(a, upper), grid, cumulative)
    return numer / denom


def simulate_coverage(
    chrom_length: int,
    depth: float,
    dist: LengthDistribution,
    min_length: float,
    window: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Simulate per-base coverage in the terminal window of a linear molecule.

    Fragments are placed by drawing lengths from ``dist``, discarding those
    below ``min_length`` (size selection acts by removal, so terminal coverage
    is depleted rather than redistributed), and dropping uniform start points.
    Only the first ``window`` bases are tallied, using a difference array, so
    the cost is O(fragments + window) regardless of molecule size.

    Parameters
    ----------
    chrom_length : int
        Length of the simulated linear molecule in bases. Should be much
        larger than ``window`` and the typical fragment length so the tallied
        end behaves as semi-infinite.
    depth : float
        Target interior sequencing depth after size selection.
    dist : FixedLength or LognormalLengths
        Fragment-length distribution before size selection.
    min_length : float
        Size-selection cutoff in bases.
    window : int
        Number of terminal bases to tally (positions 0..window-1).
    rng : np.random.Generator
        Source of randomness; pass a seeded generator for reproducibility.

    Returns
    -------
    np.ndarray
        Integer coverage per position, shape ``(window,)``.

    Raises
    ------
    ValueError
        If the window does not fit in the molecule, or if ``min_length``
        discards (essentially) every fragment.
    """
    if window > chrom_length:
        raise ValueError(f'Window {window} exceeds molecule length {chrom_length}')

    trunc_mean = dist.truncated_mean(min_length)
    n_target = int(round(depth * chrom_length / trunc_mean))

    # Size selection by rejection: oversample by the inverse survival
    # probability so the post-filter fragment count lands near the target.
    if min_length > 0:
        surv = float(dist.survival(np.array([min_length]))[0])
        surv = max(surv, _MIN_SURVIVAL)
    else:
        surv = 1.0
    n_draw = int(math.ceil(n_target / surv * 1.02)) + 16
    lengths = dist.sample(rng, n_draw)
    lengths = lengths[lengths >= min_length][:n_target]
    logger.debug(
        'Simulating %d fragments (drew %d, cutoff %s bp) on a %d bp molecule '
        'for interior depth %.1f',
        lengths.size,
        n_draw,
        min_length,
        chrom_length,
        depth,
    )

    # Uniform placement of each fragment within the molecule. Fragments longer
    # than the molecule cannot be placed; clip their span to the whole molecule.
    span = np.maximum(chrom_length - lengths, 0.0)
    starts = np.floor(rng.random(lengths.size) * (span + 1.0)).astype(np.int64)
    ends = np.minimum(starts + np.floor(lengths).astype(np.int64), chrom_length)

    # Difference-array tally restricted to the terminal window: a fragment
    # contributes +1 at its start and -1 past its end, then a cumulative sum
    # turns the deltas into per-base coverage.
    diff = np.zeros(window + 1, dtype=np.int64)
    in_window = starts < window
    np.add.at(diff, starts[in_window], 1)
    np.add.at(diff, np.minimum(ends[in_window], window), -1)
    return np.cumsum(diff[:window])


@dataclass(frozen=True)
class BootstrapResult:
    """
    Bootstrap summary of simulated terminal coverage.

    Attributes
    ----------
    x : np.ndarray
        Distances from the terminus, shape ``(window,)``.
    mean : np.ndarray
        Mean simulated coverage per position across replicates.
    lo : np.ndarray
        Lower confidence bound per position.
    hi : np.ndarray
        Upper confidence bound per position.
    analytic : np.ndarray
        Analytic expected coverage ``D * c(x)`` on the same positions.
    """

    x: np.ndarray
    mean: np.ndarray
    lo: np.ndarray
    hi: np.ndarray
    analytic: np.ndarray


def bootstrap_coverage(
    chrom_length: int,
    depth: float,
    dist: LengthDistribution,
    min_length: float,
    window: int,
    n_reps: int = 200,
    seed: int = 42,
    ci: float = 0.95,
) -> BootstrapResult:
    """
    Replicate the coverage simulation and summarise a confidence band.

    Runs ``n_reps`` independent simulations (each an independent sequencing
    experiment) and reports the per-position mean and percentile band, next
    to the analytic expectation for the same parameters.

    Parameters
    ----------
    chrom_length : int
        Length of the simulated linear molecule in bases.
    depth : float
        Target interior sequencing depth after size selection.
    dist : FixedLength or LognormalLengths
        Fragment-length distribution before size selection.
    min_length : float
        Size-selection cutoff in bases.
    window : int
        Number of terminal bases to tally.
    n_reps : int, optional
        Number of independent simulation replicates.
    seed : int, optional
        Seed for the replicate generators (spawned per replicate, so results
        are reproducible and replicates are independent).
    ci : float, optional
        Central confidence level for the percentile band, e.g. 0.95.

    Returns
    -------
    BootstrapResult
        Positions, per-position mean, band bounds, and the analytic curve.
    """
    root = np.random.default_rng(seed)
    tail = 100.0 * (1.0 - ci) / 2.0
    reps = np.empty((n_reps, window), dtype=np.int64)
    for i, child in enumerate(root.spawn(n_reps)):
        reps[i] = simulate_coverage(
            chrom_length, depth, dist, min_length, window, child
        )
    logger.info(
        'Bootstrapped %d coverage replicates (depth %.1f, cutoff %s bp, %d bp window)',
        n_reps,
        depth,
        min_length,
        window,
    )

    x = np.arange(window, dtype=float)
    return BootstrapResult(
        x=x,
        mean=reps.mean(axis=0),
        lo=np.percentile(reps, tail, axis=0),
        hi=np.percentile(reps, 100.0 - tail, axis=0),
        analytic=depth * expected_relative_coverage(x, dist, min_length),
    )
