"""
Backwards-compatible alias for :mod:`teloclip.core.streaming_analysis`.

This module moved to ``teloclip.core.streaming_analysis`` when the package was reorganised
into ``core``/``io``/``report`` subpackages. The names below are re-exported so
existing imports of ``teloclip.streaming_analysis`` keep working. New code should import
from ``teloclip.core.streaming_analysis`` directly; this shim may be removed in a future
major release.
"""

from teloclip.core.streaming_analysis import (
    ExtensionResult,
    collect_contig_overhangs_streaming,
    process_single_contig_extension,
    stream_contigs_for_extension,
)

__all__ = [
    'ExtensionResult',
    'collect_contig_overhangs_streaming',
    'process_single_contig_extension',
    'stream_contigs_for_extension',
]
