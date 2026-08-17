"""
Backwards-compatible alias for the HTML report modules.

``teloclip.html_report`` was split when the package was reorganised: document
assembly moved to :mod:`teloclip.report.html`, the per-read alignment panels to
:mod:`teloclip.report.panels`, the inline SVG charts to
:mod:`teloclip.report.charts` and the stylesheet to :mod:`teloclip.report.css`.
The names below are re-exported so existing imports of ``teloclip.html_report``
keep working. New code should import from the specific module; this shim may be
removed in a future major release.
"""

from teloclip.report.html import render_html_report
from teloclip.report.panels import (
    EndPanel,
    ReadRow,
    build_end_panel,
    motif_pattern,
    render_contig_panels,
)

__all__ = [
    'EndPanel',
    'ReadRow',
    'build_end_panel',
    'motif_pattern',
    'render_contig_panels',
    'render_html_report',
]
