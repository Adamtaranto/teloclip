"""
Unit tests for the extend HTML report.

The layout arithmetic is the part most likely to be silently wrong: a
one-column error in either direction still produces a plausible-looking
alignment, so the offsets are pinned explicitly here.
"""

import re

from teloclip.analysis import ContigStats, OverhangInfo
from teloclip.html_report import build_end_panel, render_html_report

CONTIG_LEN = 1000


def make_overhang(
    read: str,
    is_left: bool,
    clip: str,
    anchor: str,
    trim: int,
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

    Returns
    -------
    OverhangInfo
        A populated record consistent with the requested geometry.
    """
    net_gain = len(clip) - trim
    if is_left:
        aln_pos, aln_end = trim + 1, CONTIG_LEN
    else:
        aln_pos, aln_end = 1, CONTIG_LEN - trim

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
        anchor_seq=anchor,
    )


def panel_for(is_left: bool, overhangs, reference='N' * 20, selected=None):
    """
    Build an EndPanel for the given overhangs.

    Parameters
    ----------
    is_left : bool
        Which contig end.
    overhangs : list
        Overhang records.
    reference : str, optional
        Terminal contig window.
    selected : str, optional
        Name of the read used for the extension.

    Returns
    -------
    EndPanel
        The prepared panel.
    """
    stats = ContigStats('contig1', CONTIG_LEN)
    return build_end_panel(
        contig='contig1',
        end='left' if is_left else 'right',
        stats=stats,
        overhangs=overhangs,
        terminal_sequence=reference,
        selected_read=selected,
        max_reads=25,
        max_overhang=300,
        flagged=False,
    )


class TestLeftEndLayout:
    """Column offsets at the 5' terminus."""

    def test_flush_read_starts_at_negative_net_gain(self):
        """A clip flush with the terminus starts net_gain columns outside it."""
        oh = make_overhang('r1', True, clip='A' * 10, anchor='C' * 5, trim=0)
        panel = panel_for(True, [oh])

        row = panel.rows[0]
        assert row.net_gain == 10
        assert row.trim == 0
        # Clip occupies offsets -10..-1, anchor starts at offset 0.
        assert row.start == -10

    def test_trimmed_read_is_shifted_inward(self):
        """A read aligning 3 bases in adds 3 fewer bases and starts 3 later."""
        oh = make_overhang('r1', True, clip='A' * 10, anchor='C' * 5, trim=3)
        panel = panel_for(True, [oh])

        row = panel.rows[0]
        assert row.net_gain == 7
        assert row.trim == 3
        # Clip spans offsets -7..2; the anchor then begins at offset 3.
        assert row.start == -7
        assert row.start + len(row.clip) == 3

    def test_reference_row_starts_at_the_terminus(self):
        """The contig reference begins at offset 0, the terminal base."""
        oh = make_overhang('r1', True, clip='A' * 10, anchor='C' * 5, trim=0)
        panel = panel_for(True, [oh], reference='G' * 20)

        assert panel.reference_start == 0

    def test_long_clip_is_truncated_at_the_far_end(self):
        """Truncation keeps the terminus-adjacent bases, and re-anchors."""
        oh = make_overhang('r1', True, clip='A' * 500, anchor='C' * 5, trim=2)
        panel = panel_for(True, [oh])

        row = panel.rows[0]
        assert len(row.clip) == 300
        # The clip must still meet the anchor at offset `trim`.
        assert row.start + len(row.clip) == row.trim


class TestRightEndLayout:
    """Column offsets at the 3' terminus, mirroring the left end."""

    def test_flush_read_anchor_ends_at_the_terminus(self):
        """A flush alignment's last anchored base sits at offset 0."""
        oh = make_overhang('r1', False, clip='A' * 10, anchor='C' * 5, trim=0)
        panel = panel_for(False, [oh])

        row = panel.rows[0]
        assert row.net_gain == 10
        # Anchor spans offsets -4..0, so the clip begins at offset 1.
        assert row.start == -4
        assert row.start + len(row.anchor) - 1 == 0

    def test_trimmed_read_is_shifted_inward(self):
        """A read ending 3 bases short has its anchor end at offset -3."""
        oh = make_overhang('r1', False, clip='A' * 10, anchor='C' * 5, trim=3)
        panel = panel_for(False, [oh])

        row = panel.rows[0]
        assert row.net_gain == 7
        assert row.start == -7
        assert row.start + len(row.anchor) - 1 == -3

    def test_reference_row_ends_at_the_terminus(self):
        """The contig reference window ends at offset 0."""
        oh = make_overhang('r1', False, clip='A' * 10, anchor='C' * 5, trim=0)
        panel = panel_for(False, [oh], reference='G' * 20)

        assert panel.reference_start == -19
        assert panel.reference_start + len(panel.reference) - 1 == 0

    def test_long_clip_is_truncated_at_the_far_end(self):
        """Truncation keeps the terminus-adjacent bases."""
        oh = make_overhang('r1', False, clip='A' * 500, anchor='C' * 5, trim=0)
        panel = panel_for(False, [oh])

        assert len(panel.rows[0].clip) == 300


class TestEndsAreMirrorImages:
    """The two ends must be laid out symmetrically."""

    def test_equivalent_reads_are_equidistant_from_the_terminus(self):
        """A read and its mirror sit the same distance outside the contig."""
        left = panel_for(
            True, [make_overhang('r', True, 'A' * 12, 'C' * 6, trim=2)]
        ).rows[0]
        right = panel_for(
            False, [make_overhang('r', False, 'A' * 12, 'C' * 6, trim=2)]
        ).rows[0]

        assert left.net_gain == right.net_gain == 10
        assert left.trim == right.trim == 2

        # Left: the clip reaches 10 columns before offset 0.
        assert left.start == -10
        # Right: the clip reaches 10 columns past offset 0.
        assert right.start + len(right.anchor) + len(right.clip) - 1 == 10


class TestPanelSelection:
    """Read ordering, capping and selection marking."""

    def test_reads_are_ranked_by_net_gain(self):
        """The read contributing most sequence is shown first."""
        overhangs = [
            make_overhang('small', True, 'A' * 5, 'C' * 5, trim=0),
            make_overhang('big', True, 'A' * 50, 'C' * 5, trim=0),
            make_overhang('mid', True, 'A' * 20, 'C' * 5, trim=0),
        ]
        panel = panel_for(True, overhangs)

        assert [r.read_name for r in panel.rows] == ['big', 'mid', 'small']

    def test_read_cap_is_applied_and_total_retained(self):
        """Only max_reads rows render, but the true total is kept for the note."""
        overhangs = [
            make_overhang(f'r{i}', True, 'A' * (i + 2), 'C' * 5, trim=0)
            for i in range(40)
        ]
        stats = ContigStats('contig1', CONTIG_LEN)
        panel = build_end_panel(
            contig='contig1',
            end='left',
            stats=stats,
            overhangs=overhangs,
            terminal_sequence='N' * 20,
            selected_read=None,
            max_reads=10,
            max_overhang=300,
            flagged=False,
        )

        assert len(panel.rows) == 10
        assert panel.total_reads == 40

    def test_selected_read_is_marked(self):
        """The read used for the extension is flagged among the candidates."""
        overhangs = [
            make_overhang('winner', True, 'A' * 50, 'C' * 5, trim=0),
            make_overhang('other', True, 'A' * 20, 'C' * 5, trim=0),
        ]
        panel = panel_for(True, overhangs, selected='winner')

        picked = [r.read_name for r in panel.rows if r.selected]
        assert picked == ['winner']

    def test_no_overhangs_yields_no_panel(self):
        """An end with no supporting reads is skipped entirely."""
        assert panel_for(True, []) is None


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
            'terminal_sequences': {'contig1': ('G' * 20, 'C' * 20)},
            'selected_reads': {'contig1': {'left': 'readL', 'right': 'readR'}},
            'total_contigs': 1,
            'dry_run': False,
            'motifs': ['TTAGGG'],
            'version': '0.0.0-test',
            'command': 'teloclip extend ...',
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
        html = self._report(motifs=[])

        assert '<mark>' not in html

    def test_read_names_are_available_on_hover(self):
        """Each read row carries a title attribute naming the read."""
        html = self._report()

        assert 'title="readL' in html
        assert 'title="readR' in html

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
            terminal_sequences={'<script>x</script>': ('G' * 20, 'C' * 20)},
            selected_reads={},
        )

        body = html[html.index('<body>') :]
        # The injected markup must appear only in escaped form. (The report's
        # own inline script is the one legitimate <script> tag, so match on the
        # injected payload rather than the tag name.)
        assert '<script>x</script>' not in body
        assert '<img src=x>' not in body
        assert '&lt;script&gt;x&lt;/script&gt;' in body
        assert '&lt;img src=x&gt;' in body

    def test_dry_run_is_stated(self):
        """A dry run says so, so the report is not mistaken for applied work."""
        html = self._report(dry_run=True)

        assert 'Dry run' in html

    def test_empty_assembly_renders(self):
        """A run with no overhangs still produces a valid document."""
        html = self._report(
            stats_dict={},
            extensions_applied={},
            terminal_sequences={},
            selected_reads={},
        )

        assert html.startswith('<!doctype html>')
        assert 'No contigs were extended' in html
