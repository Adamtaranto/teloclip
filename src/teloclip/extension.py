"""
Backwards-compatible alias for :mod:`teloclip.core.extension`.

This module moved to ``teloclip.core.extension`` when the package was reorganised
into ``core``/``io``/``report`` subpackages. The names below are re-exported so
existing imports of ``teloclip.extension`` keep working. New code should import
from ``teloclip.core.extension`` directly; this shim may be removed in a future
major release.
"""

from teloclip.core.extension import (
    apply_contig_extension,
    calculate_extension_position,
    extend_contig,
    trim_contig_end,
    validate_extension,
)

__all__ = [
    'apply_contig_extension',
    'calculate_extension_position',
    'extend_contig',
    'trim_contig_end',
    'validate_extension',
]
