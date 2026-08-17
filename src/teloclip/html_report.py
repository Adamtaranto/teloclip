"""
Backwards-compatible alias for :mod:`teloclip.report.html`.

This module moved to ``teloclip.report.html`` when the package was reorganised
into ``core``/``io``/``report`` subpackages. The names below are re-exported so
existing imports of ``teloclip.html_report`` keep working. New code should import
from ``teloclip.report.html`` directly; this shim may be removed in a future
major release.
"""

from teloclip.report.html import (
    EndPanel,
    ReadRow,
    build_end_panel,
    render_contig_panels,
    render_html_report,
)

__all__ = [
    'EndPanel',
    'ReadRow',
    'build_end_panel',
    'render_contig_panels',
    'render_html_report',
]
