"""
Unit tests for the extend HTML report.

The layout arithmetic is the part most likely to be silently wrong: a
one-column error in either direction still produces a plausible-looking
alignment, so the offsets are pinned explicitly here.
"""

import re

from teloclip.analysis import ContigStats, OverhangInfo
from teloclip.html_report import render_contig_panels, render_html_report

CONTIG_LEN = 1000


def make_overhang(
    read: str,
    is_left: bool,
    clip: str,
    anchor: str,
    trim: int,
    spans_both_ends: bool = False,
) -> OverhangInfo:
    """
    Build an OverhangInfo positioned a given distance inside the contig.

    Parameters
    ----------
    read : str
        Read name.
    is_left : bool
        Whether this is a left-end overhang.
    clip : str
        Soft-clip sequence.
    anchor : str
        Aligned sequence adjacent to the clip.
    trim : int
        Contig bases between the terminus and the alignment.
    spans_both_ends : bool, optional
        Whether the read also overhangs the other end (default: False).

    Returns
    -------
    OverhangInfo
        A populated record consistent with the requested geometry.
    """
    net_gain = len(clip) - trim
    if is_left:
        aln_pos, aln_end = trim + 1, trim + len(anchor)
        read_seq = clip + anchor
        cigar = f'{len(clip)}S{len(anchor)}M'
    else:
        aln_end = CONTIG_LEN - trim
        aln_pos = aln_end - len(anchor) + 1
        read_seq = anchor + clip
        cigar = f'{len(anchor)}M{len(clip)}S'

    return OverhangInfo(
        sequence=clip,
        length=len(clip),
        alignment_pos=aln_pos,
        alignment_end=aln_end,
        read_name=read,
        is_left=is_left,
        clip_length=len(clip),
        anchor_length=len(anchor),
        contig_name='contig1',
        net_gain=net_gain,
        read_seq=read_seq,
        cigar=cigar,
        mapq=60,
        flag=0,
        spans_both_ends=spans_both_ends,
    )


def panels_for(stats, selected=None, flagged=False, motifs=()):
    """
    Render both alignment blocks for a contig.

    Parameters
    ----------
    stats : ContigStats
        Populated contig statistics.
    selected : dict, optional
        ``{'left': read, 'right': read}`` for applied extensions.
    flagged : bool, optional
        Whether both ends are flagged for anomalous coverage.
    motifs : sequence, optional
        Motifs to highlight.

    Returns
    -------
    list of str
        Rendered HTML blocks.
    """
    return render_contig_panels(
        contig=stats.contig_name,
        stats=stats,
        terminal_sequences=('G' * 200, 'C' * 200),
        selected_reads=selected or {},
        flagged_left=flagged,
        flagged_right=flagged,
        motifs=motifs,
        max_reads=25,
        max_overhang=300,
        min_window=60,
        window_cap=200,
    )


class TestRenderedDocument:
    """End-to-end rendering."""

    def _report(self, **overrides):
        """
        Render a report over a small two-contig assembly.

        Parameters
        ----------
        **overrides
            Keyword arguments overriding the defaults.

        Returns
        -------
        str
            The rendered HTML.
        """
        stats = ContigStats('contig1', CONTIG_LEN)
        stats.left_overhangs = [
            make_overhang('readL', True, 'TTAGGG' * 4, 'ACGT' * 5, trim=2)
        ]
        stats.right_overhangs = [
            make_overhang('readR', False, 'CCCTAA' * 4, 'ACGT' * 5, trim=0)
        ]

        kwargs = {
            'stats_dict': {'contig1': stats},
            'extensions_applied': {
                'contig1': {
                    'original_length': CONTIG_LEN,
                    'final_length': CONTIG_LEN + 46,
                    'has_left_extension': True,
                    'left_net_gain': 22,
                    'left_trim_length': 2,
                    'left_read_name': 'readL',
                    'has_right_extension': True,
                    'right_net_gain': 24,
                    'right_trim_length': 0,
                    'right_read_name': 'readR',
                }
            },
            'anomalous': {'left_outliers': [], 'right_outliers': []},
            'excluded_contigs': [],
            'warnings': [],
            'panels': panels_for(
                stats,
                selected={'left': 'readL', 'right': 'readR'},
                motifs=['TTAGGG'],
            ),
            'total_contigs': 1,
            'dry_run': False,
            'motifs': ['TTAGGG'],
            'version': '0.0.0-test',
            'command': 'teloclip extend reads.bam ref.fa --html-report out.html',
            'generated': '2026-07-31 12:00:00 UTC',
        }
        kwargs.update(overrides)
        return render_html_report(**kwargs)

    def test_document_is_well_formed_and_self_contained(self):
        """The report is one HTML file with no external references."""
        html = self._report()

        assert html.startswith('<!doctype html>')
        assert html.rstrip().endswith('</html>')
        # No CDN scripts, stylesheets, fonts or images. The one script is
        # inline, so the file still opens with nothing else present.
        assert not re.search(r'(src|href)="https?://', html)
        assert '@import' not in html
        assert not re.search(r'<script[^>]*\ssrc=', html)

    def test_summary_and_sections_present(self):
        """The expected sections render."""
        html = self._report()

        for heading in (
            'Summary',
            'Extensions',
            'Overhang depth across contigs',
            'Overhang alignments',
        ):
            assert heading in html

    def test_extension_row_arithmetic_reconciles(self):
        """Original plus net equals final in the rendered table."""
        html = self._report()
        row = re.search(r'<tr><td>contig1</td>(.*?)</tr>', html, re.S).group(1)
        cells = re.findall(r'<td[^>]*>([^<]*)</td>', row)
        num = lambda s: int(s.replace(',', '').replace('+', '').replace('–', '0'))  # noqa: E731

        original, final, left, right, total = (num(cells[i]) for i in range(5))
        assert original + total == final
        assert left + right == total

    def test_motifs_are_highlighted(self):
        """Requested motifs are marked in the sequence blocks."""
        html = self._report()

        assert '<mark>TTAGGG</mark>' in html

    def test_no_highlighting_without_motifs(self):
        """No marks appear when no motifs were requested."""
        stats = ContigStats('contig1', CONTIG_LEN)
        stats.left_overhangs = [
            make_overhang('readL', True, 'TTAGGG' * 4, 'ACGT' * 5, 2)
        ]

        html = self._report(motifs=[], panels=panels_for(stats, motifs=()))

        assert '<mark>' not in html

    def test_read_rows_carry_sam_metadata_for_the_hover_table(self):
        """Each read row exposes the SAM fields the tooltip renders."""
        html = self._report()

        assert 'data-read="readL"' in html
        assert 'data-read="readR"' in html
        # The tooltip is built from these, so they have to be present.
        for attr in (
            'data-cigar=',
            'data-mapq=',
            'data-flag=',
            'data-span=',
            'data-anchor=',
            'data-clip=',
            'data-gain=',
            'data-trim=',
        ):
            assert attr in html, attr

    def test_position_ruler_is_zeroed_at_the_terminus(self):
        """The x-axis ruler labels offsets from the contig end."""
        html = self._report()

        assert 'aln-ruler' in html
        assert '<i>0</i>' in html
        # Offsets either side are signed.
        assert re.search(r'<i>[+-]\d+</i>', html)
        # Ticks position in `ch` so they track the monospace sequence grid.
        assert re.search(r'class="tick" style="left:\d+ch"', html)

    def test_read_spanning_both_ends_is_indicated(self):
        """A read hanging off both ends is called out."""
        stats = ContigStats('contig1', CONTIG_LEN)
        stats.left_overhangs = [
            make_overhang(
                'spanner', True, 'AC' * 8, 'ACGT' * 10, 0, spans_both_ends=True
            )
        ]
        stats.right_overhangs = [
            make_overhang(
                'spanner', False, 'GT' * 8, 'ACGT' * 10, 0, spans_both_ends=True
            )
        ]

        html = self._report(stats_dict={'contig1': stats}, panels=panels_for(stats))

        assert 'spans both ends' in html
        assert 'data-both="yes"' in html
        assert 'overhangs BOTH ends' in html

    def test_read_not_spanning_both_ends_is_not_flagged(self):
        """The both-ends marker is absent when it does not apply."""
        html = self._report()

        assert 'data-both="no"' in html
        assert 'spans both ends' not in html

    def test_chart_points_expose_contig_on_hover_and_click(self):
        """Chart points carry the data the tooltip and pin label read."""
        html = self._report()

        assert 'data-contig="contig1"' in html
        assert 'data-depth=' in html
        # Native title survives as a no-JavaScript fallback.
        assert re.search(r'<title>contig1 — (left|right) end', html)
        # And the click handler that pins the name.
        assert "addEventListener('click'" in html

    def test_chart_has_a_table_view(self):
        """The chart is accompanied by a table equivalent."""
        html = self._report()

        assert 'Table view' in html
        assert 'Overhang reads' in html

    def test_flagged_contig_is_reported_but_not_excluded(self):
        """Anomalous coverage is surfaced with its caveat."""
        html = self._report(
            anomalous={'left_outliers': ['contig1'], 'right_outliers': []}
        )

        assert 'anomalous overhang coverage' in html.lower()
        assert 'not been excluded' in html

    def test_html_is_escaped(self):
        """Contig and read names are escaped, not injected."""
        stats = ContigStats('<script>x</script>', CONTIG_LEN)
        stats.left_overhangs = [
            make_overhang('<img src=x>', True, 'ACGT' * 4, 'ACGT' * 4, trim=0)
        ]

        html = self._report(
            stats_dict={'<script>x</script>': stats},
            extensions_applied={},
            panels=panels_for(stats),
        )

        body = html[html.index('<body>') :]
        # The injected markup must appear only in escaped form. (The report's
        # own inline script is the one legitimate <script> tag, so match on the
        # injected payload rather than the tag name.)
        assert '<script>x</script>' not in body
        assert '<img src=x>' not in body
        assert '&lt;script&gt;x&lt;/script&gt;' in body
        assert '&lt;img src=x&gt;' in body

    def test_footer_records_provenance(self):
        """Version, timestamp and the full command close the report.

        The report is an artefact people keep; without these it cannot be
        traced back to the run that produced it.
        """
        html = self._report()
        footer = re.search(r'<footer>(.*?)</footer>', html, re.S).group(1)

        assert 'teloclip 0.0.0-test' in footer
        assert '2026-07-31 12:00:00 UTC' in footer
        # The command is preserved in full, not truncated or summarised.
        assert 'teloclip extend reads.bam ref.fa --html-report out.html' in footer
        # Labelled, so the version cannot be misread as part of the command.
        for label in ('Version', 'Generated', 'Command'):
            assert f'<dt>{label}</dt>' in footer

    def test_footer_survives_missing_provenance(self):
        """A report rendered without version or command still closes cleanly."""
        html = self._report(version='', command='', generated='')
        footer = re.search(r'<footer>(.*?)</footer>', html, re.S).group(1)

        assert 'teloclip' in footer
        assert '<dt>Command</dt>' not in footer

    def test_ruler_labels_are_lifted_off_the_sequence(self):
        """Tick labels sit above the mark, not flush with the first row.

        The label is positioned against the bottom of its tick, and the tick
        occupies only the lower part of the ruler row, which is what puts air
        between the numbers and the contig sequence below.
        """
        html = self._report()

        assert '.ruler .tick i' in html
        assert 'bottom: 100%' in html
        # The tick mark starts partway down the row rather than at its top.
        assert re.search(r'\.ruler \.tick \{[^}]*top: 1\.5em', html, re.S)

    def test_dry_run_is_stated(self):
        """A dry run says so, so the report is not mistaken for applied work."""
        html = self._report(dry_run=True)

        assert 'Dry run' in html

    def test_empty_assembly_renders(self):
        """A run with no overhangs still produces a valid document."""
        html = self._report(
            stats_dict={},
            extensions_applied={},
            panels=[],
        )

        assert html.startswith('<!doctype html>')
        assert 'No contigs were extended' in html
