"""
Backwards-compatible alias for :mod:`teloclip.core.motifs`.

This module moved to ``teloclip.core.motifs`` when the package was reorganised
into ``core``/``io``/``report`` subpackages. The names below are re-exported so
existing imports of ``teloclip.motifs`` keep working. New code should import
from ``teloclip.core.motifs`` directly; this shim may be removed in a future
major release.
"""

from teloclip.core.motifs import (
    check_sequence_for_patterns,
    construct_regex_pattern,
    count_continuous_runs,
    count_regex_patterns_in_sequence,
    format_pattern_counts,
    make_fuzzy_motif_regex,
    make_motif_regex,
)

__all__ = [
    'check_sequence_for_patterns',
    'construct_regex_pattern',
    'count_continuous_runs',
    'count_regex_patterns_in_sequence',
    'format_pattern_counts',
    'make_fuzzy_motif_regex',
    'make_motif_regex',
]
