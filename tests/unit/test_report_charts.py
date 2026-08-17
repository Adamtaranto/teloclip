"""
Unit tests for the inline SVG charts in the extend HTML report.

These check the things that are cheap to get wrong and expensive to notice: an
empty assembly producing a chart with no marks, a contig end with one read
crashing a density estimate, contig names reaching the markup unescaped, and
the mark counts drifting from the data. Appearance is not tested here — that
needs a browser and an eye — but every chart's structural contract is.
"""

from html import escape
import re

import pytest

from teloclip.core.analysis import ContigStats, OverhangInfo
from teloclip.report.charts import (
    LABEL_LINE_HEIGHT,
    MIN_READS_FOR_DENSITY,
    _density_profile,
    _flag_note,
    _place_labels,
    _quantile,
    coverage_chart,
    depth_vs_length_chart,
    overhang_length_chart,
)


def make_stats(spec, contig_length=10_000):
    """
    Build a stats mapping from per-contig overhang length lists.

    Parameters
    ----------
    spec : dict
        Mapping of contig name to a ``(left_lengths, right_lengths)`` tuple.
    contig_length : int, optional
        Length to give each contig (default: 10,000).

    Returns
    -------
    dict
        Mapping of contig name to a populated ContigStats.
    """
    stats_dict = {}
    for name, (left_lengths, right_lengths) in spec.items():
        stats = ContigStats(name, contig_length)
        stats.left_overhangs = [
            OverhangInfo('A' * n, n, 1, 100, f'l{i}', True, n, 96, name, n)
            for i, n in enumerate(left_lengths)
        ]
        stats.right_overhangs = [
            OverhangInfo(
                'T' * n,
                n,
                contig_length - 100,
                contig_length,
                f'r{i}',
                False,
                n,
                96,
                name,
                n,
            )
            for i, n in enumerate(right_lengths)
        ]
        stats_dict[name] = stats
    return stats_dict


def count_marks(markup, css_class):
    """
    Count SVG groups carrying a given leading class.

    Parameters
    ----------
    markup : str
        SVG markup to search.
    css_class : str
        Leading class name, e.g. ``'pt'`` or ``'vio'``.

    Returns
    -------
    int
        Number of matching groups.
    """
    return len(re.findall(rf'<g class="{css_class}[ "]', markup))


ALL_CHARTS = [coverage_chart, overhang_length_chart, depth_vs_length_chart]


class TestHelpers:
    """The small numeric helpers the charts are built on."""

    def test_quantile_of_single_value(self):
        """A one-element sample is its own every quantile."""
        assert _quantile([42.0], 0.25) == 42.0
        assert _quantile([42.0], 0.99) == 42.0

    def test_quantile_interpolates(self):
        """Quantiles falling between observations are interpolated."""
        assert _quantile([0.0, 10.0], 0.5) == 5.0
        assert _quantile([0.0, 10.0, 20.0, 30.0], 0.5) == 15.0

    @pytest.mark.parametrize('q,expected', [(0.0, 1.0), (1.0, 4.0)])
    def test_quantile_endpoints(self, q, expected):
        """
        The extreme quantiles return the extreme observations.

        Parameters
        ----------
        q : float
            Quantile requested.
        expected : float
            Value expected back.
        """
        assert _quantile([1.0, 2.0, 3.0, 4.0], q) == expected

    def test_density_profile_is_normalised(self):
        """The profile peaks at 1.0 so violins are comparable in width."""
        profile = _density_profile([10.0, 20.0, 30.0, 40.0], 0.0, 100.0)

        assert max(profile) == pytest.approx(1.0)
        assert all(0.0 <= weight <= 1.0 for weight in profile)

    def test_density_profile_of_nothing_is_flat_zero(self):
        """An empty sample has no density anywhere, and does not divide by zero."""
        assert _density_profile([], 0.0, 100.0) == [0.0] * 14

    def test_density_profile_survives_a_zero_width_scale(self):
        """
        Every read the same length must not divide by a zero span.

        This happens in practice whenever a contig end is supported by reads
        clipped at exactly the same point.
        """
        profile = _density_profile([50.0, 50.0, 50.0], 50.0, 50.0)

        assert max(profile) == pytest.approx(1.0)

    def test_density_profile_puts_weight_near_the_data(self):
        """Mass lands in the half of the scale the observations occupy."""
        profile = _density_profile([90.0] * 10, 0.0, 100.0)

        assert sum(profile[7:]) > sum(profile[:7])

    @pytest.mark.parametrize(
        'depth,length,expected',
        [
            (True, True, 'anomalous depth and length'),
            (True, False, 'anomalous depth'),
            (False, True, 'anomalous length'),
            (False, False, '–'),
        ],
    )
    def test_flag_note_names_the_combination(self, depth, length, expected):
        """
        Each combination of the two flags is described distinctly.

        Parameters
        ----------
        depth : bool
            Whether depth was flagged.
        length : bool
            Whether length was flagged.
        expected : str
            Expected description.
        """
        assert _flag_note(depth, length) == expected


class TestEmptyInput:
    """No data means no chart, rather than an empty pair of axes."""

    @pytest.mark.parametrize('chart', ALL_CHARTS, ids=lambda f: f.__name__)
    def test_empty_stats_dict(self, chart):
        """
        An assembly with no contigs produces nothing at all.

        Parameters
        ----------
        chart : callable
            Chart function under test.
        """
        assert chart({}, [], []) == ('', '')

    @pytest.mark.parametrize(
        'chart',
        [overhang_length_chart, depth_vs_length_chart],
        ids=lambda f: f.__name__,
    )
    def test_contigs_with_no_overhangs(self, chart):
        """
        Contigs present but with no overhangs produce nothing.

        The two length-based charts have nothing to say about an end with no
        reads, so they omit the section rather than draw an empty axis. The
        depth chart still plots those ends, at zero.

        Parameters
        ----------
        chart : callable
            Chart function under test.
        """
        stats = make_stats({'contig1': ([], []), 'contig2': ([], [])})

        assert chart(stats, [], []) == ('', '')

    def test_coverage_chart_still_plots_empty_ends(self):
        """Zero depth is a value the depth chart should show, not omit."""
        stats = make_stats({'contig1': ([], []), 'contig2': ([100], [])})

        chart, table = coverage_chart(stats, [], [])

        assert chart
        # Two ends per contig, whether or not they carry reads.
        assert count_marks(chart, 'pt') == 4


class TestOverhangLengthChart:
    """The split-violin distribution chart."""

    def test_one_group_per_non_empty_end(self):
        """Ends with no reads get no shape; ends with reads get exactly one."""
        stats = make_stats(
            {
                'contig1': ([10, 20, 30, 40], [15, 25, 35, 45]),
                'contig2': ([10, 20, 30, 40], []),
            }
        )

        chart, _ = overhang_length_chart(stats, [], [])

        assert count_marks(chart, 'vio') == 3

    def test_sparse_end_falls_back_to_points(self):
        """
        Fewer reads than the density threshold are drawn individually.

        A density curve through two points describes the estimator rather than
        the data, so those ends show the reads themselves.
        """
        sparse = [10, 20]
        assert len(sparse) < MIN_READS_FOR_DENSITY
        stats = make_stats({'contig1': (sparse, [])})

        chart, _ = overhang_length_chart(stats, [], [])

        assert '<polygon class="vio-body"' not in chart
        assert chart.count('<circle class="vio-dot"') == len(sparse)

    def test_dense_end_draws_a_density_shape(self):
        """At or above the threshold the end is drawn as a violin."""
        lengths = [10, 20, 30, 40, 50, 60]
        assert len(lengths) >= MIN_READS_FOR_DENSITY
        stats = make_stats({'contig1': (lengths, [])})

        chart, _ = overhang_length_chart(stats, [], [])

        assert chart.count('<polygon class="vio-body"') == 1
        assert '<circle class="vio-dot"' not in chart

    def test_single_read_end_does_not_crash(self):
        """One read is a degenerate distribution but still has to render."""
        stats = make_stats({'contig1': ([250], [])})

        chart, table = overhang_length_chart(stats, [], [])

        assert chart
        assert '250' in table

    def test_every_end_gets_a_median_tick(self):
        """The median is marked whether the end is drawn as a shape or points."""
        stats = make_stats({'contig1': ([10, 20], [10, 20, 30, 40, 50, 60])})

        chart, _ = overhang_length_chart(stats, [], [])

        assert chart.count('class="vio-median"') == 2

    def test_flagged_end_is_marked_and_labelled(self):
        """A flagged contig carries the flag class and a direct label."""
        stats = make_stats({f'contig{i}': ([10, 20, 30, 40], []) for i in range(3)})

        chart, table = overhang_length_chart(stats, ['contig1'], [])

        assert 'vio-flag' in chart
        assert 'class="pt-label"' in chart
        assert 'is-flagged' in table

    def test_scale_clipping_is_disclosed(self):
        """
        When the axis is clipped the caption says so, and names the true maximum.

        A silently clipped axis would understate the longest overhang in the
        assembly, which is exactly the value a reader is looking for.
        """
        spec = {f'contig{i}': ([100, 105, 110, 115], []) for i in range(12)}
        spec['huge'] = ([100_000, 100_001, 100_002, 100_003], [])
        stats = make_stats(spec, contig_length=200_000)

        chart, _ = overhang_length_chart(stats, [], [])

        assert 'clipped at' in chart
        assert '100,003' in chart

    def test_no_clipping_note_when_nothing_is_clipped(self):
        """A caption that always claims clipping would be noise."""
        stats = make_stats({'contig1': ([10, 11, 12, 13], [])})

        chart, _ = overhang_length_chart(stats, [], [])

        assert 'clipped at' not in chart

    def test_table_row_per_non_empty_end(self):
        """The table twin covers every end the chart draws."""
        stats = make_stats(
            {'contig1': ([10, 20, 30, 40], [50, 60, 70, 80]), 'contig2': ([90], [])}
        )

        _, table = overhang_length_chart(stats, [], [])

        assert table.split('<tbody>')[1].count('<tr') == 3


class TestDepthVsLengthChart:
    """The linked scatter of depth against median length."""

    def test_one_point_per_non_empty_end(self):
        """Each end with reads is one point; empty ends are absent."""
        stats = make_stats(
            {
                'contig1': ([10, 20], [30, 40]),
                'contig2': ([50], []),
            }
        )

        chart, _ = depth_vs_length_chart(stats, [], [])

        assert count_marks(chart, 'pt') == 3

    def test_points_carry_the_metadata_the_tooltip_shows(self):
        """
        Every field the tooltip reads is present on the mark.

        The tooltip is built from data attributes, so a missing one degrades
        silently to a shorter tooltip rather than an error.
        """
        stats = make_stats({'contig1': ([10, 20, 30], [])}, contig_length=5_000)

        chart, _ = depth_vs_length_chart(stats, [], [])

        for attribute in (
            'data-contig',
            'data-end',
            'data-depth',
            'data-contiglength',
            'data-median',
            'data-longest',
            'data-gain',
            'data-flagged',
            'data-note',
        ):
            assert attribute in chart, f'{attribute} missing from the scatter marks'

    def test_both_ends_share_a_contig_key_for_linking(self):
        """
        Both ends of a contig carry the same data-contig value.

        This is what lets one click highlight the contig across all three
        charts; if the values diverged the linking would silently half-work.
        """
        stats = make_stats({'contig1': ([10, 20], [30, 40])})

        chart, _ = depth_vs_length_chart(stats, [], [])

        assert chart.count('data-contig="contig1"') == 2

    def test_length_flags_are_honoured_separately_from_depth(self):
        """
        An end flagged only on length is still marked as anomalous.

        The scatter is where the two measures are compared, so it has to
        reflect both sets of flags rather than only the depth ones.
        """
        stats = make_stats({'contig1': ([10, 20], []), 'contig2': ([30, 40], [])})

        chart, table = depth_vs_length_chart(
            stats, [], [], flagged_left_length=['contig2']
        )

        assert 'anomalous length' in chart
        assert 'anomalous length' in table
        assert chart.count('class="ring"') == 1

    def test_end_flagged_on_both_says_so(self):
        """The combination is named, since it is the diagnostic case."""
        stats = make_stats({'contig1': ([10, 20], [])})

        chart, _ = depth_vs_length_chart(
            stats, ['contig1'], [], flagged_left_length=['contig1']
        )

        assert 'anomalous depth and length' in chart

    def test_median_is_used_not_mean(self):
        """
        The x position summarises with the median.

        One read running far past the end would otherwise move the whole point,
        which is the failure the median exists to avoid here.
        """
        stats = make_stats({'contig1': ([10, 10, 10, 10_000], [])})

        chart, _ = depth_vs_length_chart(stats, [], [])

        assert 'data-median="10"' in chart

    def test_table_reports_best_gain(self):
        """The most a contig could gain at that end is in the table."""
        stats = make_stats({'contig1': ([10, 40, 25], [])})

        _, table = depth_vs_length_chart(stats, [], [])

        assert '40' in table


class TestEscaping:
    """Contig names come from the reference and are never trusted."""

    HOSTILE = 'contig<script>&"1'

    @pytest.mark.parametrize('chart', ALL_CHARTS, ids=lambda f: f.__name__)
    def test_contig_names_are_escaped(self, chart):
        """
        A contig name containing markup never reaches the output raw.

        FASTA headers are arbitrary text, so this is reachable with a
        deliberately named reference.

        Parameters
        ----------
        chart : callable
            Chart function under test.
        """
        stats = make_stats({self.HOSTILE: ([10, 20, 30, 40], [50, 60, 70, 80])})

        markup, table = chart(stats, [self.HOSTILE], [self.HOSTILE])

        assert '<script>' not in markup
        assert '<script>' not in table
        assert escape(self.HOSTILE) in markup


class TestSharedContract:
    """Properties every chart in the report holds to."""

    @pytest.mark.parametrize('chart', ALL_CHARTS, ids=lambda f: f.__name__)
    def test_chart_ships_a_table_twin(self, chart):
        """
        Every chart has an equivalent table, so no value is SVG-only.

        Parameters
        ----------
        chart : callable
            Chart function under test.
        """
        stats = make_stats({'contig1': ([10, 20, 30, 40], [50, 60, 70, 80])})

        markup, table = chart(stats, [], [])

        assert 'class="table-view"' in table
        assert '<table>' in table

    @pytest.mark.parametrize('chart', ALL_CHARTS, ids=lambda f: f.__name__)
    def test_chart_is_labelled_and_scrollable(self, chart):
        """
        Each chart is announced to assistive technology and scrolls in place.

        The scroll container is what keeps a wide plot from making the whole
        page scroll sideways.

        Parameters
        ----------
        chart : callable
            Chart function under test.
        """
        stats = make_stats({'contig1': ([10, 20, 30, 40], [50, 60, 70, 80])})

        markup, _ = chart(stats, [], [])

        assert 'class="chart-scroll"' in markup
        assert 'role="img"' in markup
        assert 'aria-label="' in markup

    @pytest.mark.parametrize('chart', ALL_CHARTS, ids=lambda f: f.__name__)
    def test_marks_are_focusable_and_titled(self, chart):
        """
        Marks are reachable by keyboard and carry a native title.

        The title is the no-JavaScript fallback for the tooltip.

        Parameters
        ----------
        chart : callable
            Chart function under test.
        """
        stats = make_stats({'contig1': ([10, 20, 30, 40], [50, 60, 70, 80])})

        markup, _ = chart(stats, [], [])

        assert 'tabindex="0"' in markup
        assert '<title>' in markup

    @pytest.mark.parametrize('chart', ALL_CHARTS, ids=lambda f: f.__name__)
    def test_a_large_assembly_renders(self, chart):
        """
        Two hundred contigs produce a chart rather than an exception.

        Parameters
        ----------
        chart : callable
            Chart function under test.
        """
        stats = make_stats(
            {
                f'contig{i:03d}': ([10 + i, 20 + i, 30 + i, 40 + i], [])
                for i in range(200)
            }
        )

        markup, table = chart(stats, [], [])

        assert markup
        assert table


class TestFlaggedLabelling:
    """
    Direct labels are per contig, not per contig end.

    A contig flagged at both ends is the common case — a high-copy circular
    element is deep at both — and labelling each end drew the same name twice
    a few pixels apart, where the two copies overlapped each other and the
    marks they belonged to.
    """

    def spec(self):
        """
        Build an assembly with one contig flagged at both ends.

        Returns
        -------
        dict
            Contig name to ContigStats.
        """
        spec = {f'contig{i}': ([10, 20, 30, 40], [10, 20, 30, 40]) for i in range(9)}
        spec['plasmid'] = ([50] * 80, [50] * 80)
        return make_stats(spec)

    def test_coverage_chart_labels_a_contig_once(self):
        """Both ends flagged yields one label, not two."""
        chart, _ = coverage_chart(self.spec(), ['plasmid'], ['plasmid'])

        assert chart.count('>plasmid</text>') == 1

    def test_scatter_labels_a_contig_once(self):
        """Same for the scatter, where both ends are separate points."""
        chart, _ = depth_vs_length_chart(self.spec(), ['plasmid'], ['plasmid'])

        assert chart.count('>plasmid</text>') == 1

    def test_distinct_flagged_contigs_are_each_labelled(self):
        """De-duplicating by contig must not collapse different contigs."""
        chart, _ = coverage_chart(self.spec(), ['plasmid', 'contig3'], [])

        assert chart.count('>plasmid</text>') == 1
        assert chart.count('>contig3</text>') == 1

    def test_both_ends_still_marked_even_though_labelled_once(self):
        """
        Only the label is de-duplicated; both marks keep their flag ring.

        The label is a convenience; the ring is the actual signal, and
        suppressing it on one end would understate the problem.
        """
        chart, _ = coverage_chart(self.spec(), ['plasmid'], ['plasmid'])

        assert chart.count('class="ring"') == 2


class TestLabelPlacement:
    """
    Direct labels are stacked so they never overlap.

    Contig names are long and flagged contigs cluster together, so two labels
    placed at their marks routinely collided and rendered as one unreadable
    run of characters.
    """

    def test_isolated_labels_do_not_move(self):
        """A label with nothing near it stays exactly where it was put."""
        entries = [(100.0, 50.0, 'a'), (800.0, 50.0, 'b')]

        assert _place_labels(entries) == [(100.0, 50.0, 'a'), (800.0, 50.0, 'b')]

    def test_overlapping_labels_are_separated(self):
        """Two labels at the same spot end up on different lines."""
        entries = [(100.0, 50.0, 'chr7_rdna_array'), (110.0, 50.0, 'chr8_telomere')]

        placed = _place_labels(entries)
        tops = sorted(y for _, y, _ in placed)

        assert len(set(tops)) == 2
        assert tops[1] - tops[0] >= LABEL_LINE_HEIGHT

    def test_labels_far_apart_vertically_do_not_move(self):
        """Sharing an x is fine when the labels are already on different rows."""
        entries = [(100.0, 50.0, 'alpha'), (100.0, 200.0, 'beta')]

        assert _place_labels(entries) == entries

    def test_three_crowded_labels_all_separate(self):
        """A run of collisions stacks rather than piling up on one line."""
        entries = [
            (100.0, 50.0, 'aaaaaaaaaa'),
            (105.0, 50.0, 'bbbbbbbbbb'),
            (110.0, 50.0, 'cccccccccc'),
        ]

        tops = [y for _, y, _ in _place_labels(entries)]

        assert len(set(tops)) == 3

    def test_labels_only_move_upward(self):
        """
        Labels move up, never down.

        Down would take them over the mark they belong to, which is the one
        place they must stay clear of.
        """
        entries = [(100.0, 50.0, 'alpha'), (105.0, 50.0, 'beta')]

        assert all(y <= 50.0 for _, y, _ in _place_labels(entries))

    def test_crowded_flagged_contigs_render_without_collision(self):
        """
        End to end: adjacent flagged contigs with long names get distinct rows.

        This is the case that produced an unreadable overlap in a real report.
        """
        spec = {f'contig{i}': ([10, 20, 30, 40], []) for i in range(9)}
        spec['chr7_rdna_array_collapsed'] = ([3000] * 6, [])
        spec['chr8_long_telomere_unfinished'] = ([3100] * 6, [])
        stats = make_stats(spec, contig_length=100_000)

        chart, _ = overhang_length_chart(
            stats,
            ['chr7_rdna_array_collapsed', 'chr8_long_telomere_unfinished'],
            [],
        )

        ys = re.findall(r'<text class="pt-label" x="[\d.]+" y="([-\d.]+)"', chart)
        assert len(ys) == 2
        assert len(set(ys)) == 2

    def test_stacked_labels_stay_inside_the_violin_chart(self):
        """
        The top padding leaves room for the rows the stacking can produce.

        A label pushed above the viewBox is cropped by the card edge, which is
        worse than the overlap it was avoiding.
        """
        spec = {f'contig{i}': ([10, 20, 30, 40], []) for i in range(9)}
        spec['chr7_rdna_array_collapsed'] = ([3000] * 6, [])
        spec['chr8_long_telomere_unfinished'] = ([3100] * 6, [])
        stats = make_stats(spec, contig_length=100_000)

        chart, _ = overhang_length_chart(
            stats,
            ['chr7_rdna_array_collapsed', 'chr8_long_telomere_unfinished'],
            [],
        )

        ys = [
            float(y)
            for y in re.findall(r'class="pt-label" x="[\d.]+" y="([-\d.]+)"', chart)
        ]
        assert all(y > 0 for y in ys), f'label above the top of the chart: {ys}'
