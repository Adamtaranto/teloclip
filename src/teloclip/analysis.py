"""
Backwards-compatible alias for :mod:`teloclip.core.analysis`.

This module moved to ``teloclip.core.analysis`` when the package was reorganised
into ``core``/``io``/``report`` subpackages. The names below are re-exported so
existing imports of ``teloclip.analysis`` keep working. New code should import
from ``teloclip.core.analysis`` directly; this shim may be removed in a future
major release.
"""

from teloclip.core.analysis import (
    GAIN_MARGIN_BASES,
    GAIN_MARGIN_FRACTION,
    MIN_CONTIGS_FOR_ANOMALY_FLAGGING,
    ContigStats,
    OverhangInfo,
    calculate_overhang_statistics,
    detect_homopolymer_runs,
    flag_anomalous_overhang_coverage,
    rank_overhangs_by_gain,
    select_best_overhang,
)

__all__ = [
    'ContigStats',
    'GAIN_MARGIN_BASES',
    'GAIN_MARGIN_FRACTION',
    'MIN_CONTIGS_FOR_ANOMALY_FLAGGING',
    'OverhangInfo',
    'calculate_overhang_statistics',
    'detect_homopolymer_runs',
    'flag_anomalous_overhang_coverage',
    'rank_overhangs_by_gain',
    'select_best_overhang',
]
