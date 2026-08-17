"""
Backwards-compatible alias for :mod:`teloclip.io.sam`.

This module moved to ``teloclip.io.sam`` when the package was reorganised
into ``core``/``io``/``report`` subpackages. The names below are re-exported so
existing imports of ``teloclip.samops`` keep working. New code should import
from ``teloclip.io.sam`` directly; this shim may be removed in a future
major release.
"""

from teloclip.io.sam import (
    CIGARinfo,
    EnhancedStreamingSamFilter,
    SAMinfo,
    calculate_aligned_bases,
    checkClips,
    enhanced_streaming_split_by_contig,
    lenCIGAR,
    processSamlines,
    splitCIGAR,
    validate_min_anchor,
)

__all__ = [
    'CIGARinfo',
    'EnhancedStreamingSamFilter',
    'SAMinfo',
    'calculate_aligned_bases',
    'checkClips',
    'enhanced_streaming_split_by_contig',
    'lenCIGAR',
    'processSamlines',
    'splitCIGAR',
    'validate_min_anchor',
]
