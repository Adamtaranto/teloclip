"""
Backwards-compatible alias for :mod:`teloclip.io.streaming`.

This module moved to ``teloclip.io.streaming`` when the package was reorganised
into ``core``/``io``/``report`` subpackages. The names below are re-exported so
existing imports of ``teloclip.streaming_io`` keep working. New code should import
from ``teloclip.io.streaming`` directly; this shim may be removed in a future
major release.
"""

from teloclip.io.streaming import (
    BufferedContigWriter,
    StreamingGenomeProcessor,
    validate_indexed_files,
)

__all__ = [
    'BufferedContigWriter',
    'StreamingGenomeProcessor',
    'validate_indexed_files',
]
