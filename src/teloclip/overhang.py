"""
Backwards-compatible alias for :mod:`teloclip.core.overhang`.

This module moved to ``teloclip.core.overhang`` when the package was reorganised
into ``core``/``io``/``report`` subpackages. The names below are re-exported so
existing imports of ``teloclip.overhang`` keep working. New code should import
from ``teloclip.core.overhang`` directly; this shim may be removed in a future
major release.
"""

from teloclip.core.overhang import (
    ANCHOR_OPS,
    REASON_MAX_BREAK,
    REASON_MIN_ANCHOR,
    REASON_MIN_CLIP,
    REASON_NO_CLIP,
    REF_CONSUMING,
    AlignmentEnds,
    CigarOps,
    OverhangCall,
    alignment_end,
    anchor_length,
    cigar_ops_from_string,
    cigar_ops_from_tuples,
    classify,
    classify_end,
    clip_lengths,
    ends_from_aligned_segment,
    ends_from_sam_fields,
    overhang_info_from_call,
    reference_span,
)

__all__ = [
    'ANCHOR_OPS',
    'AlignmentEnds',
    'CigarOps',
    'OverhangCall',
    'REASON_MAX_BREAK',
    'REASON_MIN_ANCHOR',
    'REASON_MIN_CLIP',
    'REASON_NO_CLIP',
    'REF_CONSUMING',
    'alignment_end',
    'anchor_length',
    'cigar_ops_from_string',
    'cigar_ops_from_tuples',
    'classify',
    'classify_end',
    'clip_lengths',
    'ends_from_aligned_segment',
    'ends_from_sam_fields',
    'overhang_info_from_call',
    'reference_span',
]
