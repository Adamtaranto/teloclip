"""
Backwards-compatible alias for :mod:`teloclip.report.layout`.

This module moved to ``teloclip.report.layout`` when the package was reorganised
into ``core``/``io``/``report`` subpackages. The names below are re-exported so
existing imports of ``teloclip.alignment_layout`` keep working. New code should import
from ``teloclip.report.layout`` directly; this shim may be removed in a future
major release.
"""

from teloclip.report.layout import (
    BLANK,
    GAP,
    Column,
    Placed,
    build_columns,
    place_read,
    render_reference,
    render_row,
    terminus_column,
)

__all__ = [
    'BLANK',
    'Column',
    'GAP',
    'Placed',
    'build_columns',
    'place_read',
    'render_reference',
    'render_row',
    'terminus_column',
]
