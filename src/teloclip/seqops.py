"""
Backwards-compatible alias for :mod:`teloclip.core.seqops`.

This module moved to ``teloclip.core.seqops`` when the package was reorganised
into ``core``/``io``/``report`` subpackages. The names below are re-exported so
existing imports of ``teloclip.seqops`` keep working. New code should import
from ``teloclip.core.seqops`` directly; this shim may be removed in a future
major release.
"""

from teloclip.core.seqops import (
    addRevComplement,
    filterList,
    isMotifInClip,
    makeMask,
    read_fai,
    revComp,
)

__all__ = [
    'addRevComplement',
    'filterList',
    'isMotifInClip',
    'makeMask',
    'read_fai',
    'revComp',
]
