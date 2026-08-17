"""
Backwards-compatible alias for :mod:`teloclip.report.text`.

This module moved to ``teloclip.report.text`` when the package was reorganised
into ``core``/``io``/``report`` subpackages. The names below are re-exported so
existing imports of ``teloclip.reporting`` keep working. New code should import
from ``teloclip.report.text`` directly; this shim may be removed in a future
major release.
"""

from teloclip.report.text import (
    fmt_delta,
    fmt_delta_html,
    fmt_float,
    fmt_int,
    histogram,
    kv_table,
    md_table,
)

__all__ = [
    'fmt_delta',
    'fmt_delta_html',
    'fmt_float',
    'fmt_int',
    'histogram',
    'kv_table',
    'md_table',
]
