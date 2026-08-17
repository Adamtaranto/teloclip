"""
Statistical model of coverage decay towards linear contig ends.

This subpackage is optional: it needs numpy (and matplotlib for the figure
helpers), which are not runtime dependencies of teloclip. Install them with
the ``model`` extra.
"""

try:
    import numpy  # noqa: F401
except ImportError as exc:  # pragma: no cover - exercised only without numpy
    raise ImportError(
        "teloclip.model requires the 'model' extra: pip install 'teloclip[model]'"
    ) from exc

from teloclip.model.decay import (
    BootstrapResult,
    FixedLength,
    LognormalLengths,
    bootstrap_coverage,
    expected_relative_coverage,
    expected_relative_coverage_quadrature,
    simulate_coverage,
)

__all__ = [
    'BootstrapResult',
    'FixedLength',
    'LognormalLengths',
    'bootstrap_coverage',
    'expected_relative_coverage',
    'expected_relative_coverage_quadrature',
    'simulate_coverage',
]
