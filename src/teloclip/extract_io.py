"""
Backwards-compatible alias for :mod:`teloclip.io.extract`.

This module moved to ``teloclip.io.extract`` when the package was reorganised
into ``core``/``io``/``report`` subpackages. The names below are re-exported so
existing imports of ``teloclip.extract_io`` keep working. New code should import
from ``teloclip.io.extract`` directly; this shim may be removed in a future
major release.
"""

from teloclip.io.extract import (
    BufferedContigWriter,
    EfficientSequenceWriter,
    ExtractionStats,
    MultiFileSequenceWriter,
    StreamingGenomeProcessor,
    create_fasta_index,
    write_sequences_streaming,
)

__all__ = [
    'BufferedContigWriter',
    'EfficientSequenceWriter',
    'ExtractionStats',
    'MultiFileSequenceWriter',
    'StreamingGenomeProcessor',
    'create_fasta_index',
    'write_sequences_streaming',
]
