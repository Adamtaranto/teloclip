"""Unit tests for the `teloclip extend` statistics report."""

import pytest

from teloclip.analysis import ContigStats, OverhangInfo
from teloclip.commands.extend import generate_extension_report


def make_overhang(length, is_left, contig='contig01', read='read1'):
    """
    Build a minimal OverhangInfo for report tests.

    Parameters
    ----------
    length : int
        Overhang length in bases.
    is_left : bool
        Whether the overhang is at the left end of the contig.
    contig : str, optional
        Contig name.
    read : str, optional
        Read name.

    Returns
    -------
    OverhangInfo
        A populated overhang record.
    """
    return OverhangInfo(
        sequence='A' * length,
        length=length,
        alignment_pos=0,
        alignment_end=length,
        read_name=read,
        is_left=is_left,
        clip_length=length,
        anchor_length=100,
        contig_name=contig,
    )


def cells(line):
    """
    Split a rendered Markdown table row into stripped cell values.

    Parameters
    ----------
    line : str
        A single rendered table line.

    Returns
    -------
    list of str
        The cell contents, with padding removed.
    """
    return [cell.strip() for cell in line.strip().strip('|').split('|')]


def find_row(report, first_cell):
    """
    Find the first table row whose leading cell matches.

    Parameters
    ----------
    report : str
        The rendered report.
    first_cell : str
        Expected value of the row's first cell.

    Returns
    -------
    list of str
        The matching row's stripped cells.

    Raises
    ------
    AssertionError
        If no row matches.
    """
    for line in report.split('\n'):
        if line.startswith('|'):
            row = cells(line)
            if row and row[0] == first_cell:
                return row
    raise AssertionError(f'no row starting with {first_cell!r}')


@pytest.fixture
def stats_dict():
    """
    Provide overhang statistics for a single contig with both ends supported.

    Returns
    -------
    dict
        Mapping of contig name to ContigStats.
    """
    return {
        'contig01': ContigStats(
            contig_name='contig01',
            contig_length=1000,
            left_overhangs=[
                make_overhang(50, True, read='readL1'),
                make_overhang(60, True, read='readL2'),
            ],
            right_overhangs=[make_overhang(80, False, read='readR1')],
        )
    }


@pytest.fixture
def both_ends_extension():
    """
    Provide an extension record for a contig extended at both ends.

    Returns
    -------
    dict
        Mapping of contig name to extension info.
    """
    return {
        'contig01': {
            'original_length': 1000,
            'final_length': 1135,
            'has_left_extension': True,
            'left_overhang_length': 60,
            'left_read_name': 'readL2',
            'left_trim_length': 5,
            'has_right_extension': True,
            'right_overhang_length': 80,
            'right_read_name': 'readR1',
            'right_trim_length': 0,
        }
    }


@pytest.fixture
def empty_outliers():
    """
    Provide an empty outlier mapping.

    Returns
    -------
    dict
        Mapping with empty left and right outlier lists.
    """
    return {'left_outliers': [], 'right_outliers': []}


def test_report_has_expected_sections(stats_dict, both_ends_extension, empty_outliers):
    """The report opens with a title and summary, then the detail tables."""
    report = generate_extension_report(
        stats_dict,
        both_ends_extension,
        empty_outliers,
        {},
        [],
        [],
    )
    assert report.startswith('# Teloclip Extend Report')
    assert '## Summary' in report
    assert '## Extensions Applied' in report
    assert '## Per-Contig Overhang Support' in report


def test_summary_counts_both_ends(stats_dict, both_ends_extension, empty_outliers):
    """A contig extended at both ends contributes two ends and the summed bases."""
    report = generate_extension_report(
        stats_dict,
        both_ends_extension,
        empty_outliers,
        {},
        [],
        [],
        total_contigs=1,
    )
    assert 'Contig ends extended' in report
    assert '2 of 2' in report
    assert '| Net bases gained' in report
    # 60 grafted less 5 trimmed on the left, plus 80 on the right.
    assert '135' in report
    assert '| Bases trimmed back' in report
    assert '| Raw overhang bases grafted' in report
    assert '140' in report  # 60 + 80 before trimming


def test_extensions_table_reports_both_ends(
    stats_dict, both_ends_extension, empty_outliers
):
    """Both left and right extension lengths appear, not just one."""
    report = generate_extension_report(
        stats_dict,
        both_ends_extension,
        empty_outliers,
        {},
        [],
        [],
    )
    row = find_row(report, 'contig01')
    assert row[1] == '1,000'
    assert row[2] == '1,135'
    # Net of the 5 bases trimmed at the left end, so the row reconciles:
    # 1,000 + 135 = 1,135.
    assert row[3] == '+55'
    assert row[4] == '+80'
    assert row[5] == '+135'
    assert row[6] == 'readL2'
    assert row[7] == 'readR1'
    assert row[8] == '5'


def test_extension_row_arithmetic_reconciles(
    stats_dict, both_ends_extension, empty_outliers
):
    """Original bp plus the net columns equals Final bp.

    Reporting raw clip lengths instead of net gain made this identity false
    wherever an end needed trimming, so the report could not be checked
    against the sequences it described.
    """
    report = generate_extension_report(
        stats_dict,
        both_ends_extension,
        empty_outliers,
        {},
        [],
        [],
    )
    row = find_row(report, 'contig01')

    def as_int(cell):
        """
        Parse a formatted report cell into an integer.

        Parameters
        ----------
        cell : str
            Cell text, possibly carrying thousands separators and a sign.

        Returns
        -------
        int
            The numeric value.
        """
        return int(cell.replace(',', '').replace('+', ''))

    assert as_int(row[1]) + as_int(row[5]) == as_int(row[2])
    assert as_int(row[3]) + as_int(row[4]) == as_int(row[5])


def test_motif_gain_is_attributed_per_end(
    stats_dict, both_ends_extension, empty_outliers
):
    """
    Motif gain is computed per contig end.

    Regression test: the previous report guessed a single end from the legacy
    `is_left` flag, charging a both-ends contig's entire gain to the left end.
    """
    terminal = {'contig01': {'left': {'TTAGGG': 1}, 'right': {'TTAGGG': 2}}}
    post = {'contig01': {'left': {'TTAGGG': 6}, 'right': {'TTAGGG': 12}}}

    report = generate_extension_report(
        stats_dict,
        both_ends_extension,
        empty_outliers,
        {},
        [],
        [],
        terminal_motif_counts=terminal,
        post_motif_counts=post,
        screen_terminal_bases=1000,
    )

    assert '## Telomere Motif Analysis' in report
    motif_rows = [
        cells(line)
        for line in report.split('\n')
        if line.startswith('| contig01 |') and ' TTAGGG ' in line
    ]
    left_row = next(r for r in motif_rows if r[1] == 'left')
    right_row = next(r for r in motif_rows if r[1] == 'right')
    # Windows are the screening window plus that end's net extension, which is
    # how much longer the sequence actually is there.
    assert left_row[3] == '1,055'
    assert left_row[4:] == ['1', '6', '+5']
    assert right_row[3] == '1,080'
    assert right_row[4:] == ['2', '12', '+10']
    assert 'Total motif gain: +15' in report


def test_motif_section_omitted_without_counts(
    stats_dict, both_ends_extension, empty_outliers
):
    """No motif section is emitted when --count-motifs was not used."""
    report = generate_extension_report(
        stats_dict,
        both_ends_extension,
        empty_outliers,
        {},
        [],
        [],
    )
    assert '## Telomere Motif Analysis' not in report
    assert 'Total motif gain' not in report


def test_dry_run_title_and_wording(stats_dict, both_ends_extension, empty_outliers):
    """Dry runs are labelled and do not carry the polishing note."""
    report = generate_extension_report(
        stats_dict,
        both_ends_extension,
        empty_outliers,
        {},
        [],
        [],
        dry_run=True,
    )
    assert report.startswith('# Teloclip Extend Report (dry run)')
    assert '## Extensions That Would Be Applied' in report
    assert 'should be polished' not in report


def test_empty_sections_are_omitted(stats_dict, empty_outliers):
    """Sections with no content are skipped rather than left as bare headings."""
    report = generate_extension_report(
        stats_dict,
        {},
        empty_outliers,
        {},
        [],
        [],
    )
    assert '## Extensions Applied' not in report
    assert '## Excluded Contigs' not in report
    assert '## Warnings' not in report


def test_excluded_and_warnings_rendered(stats_dict, empty_outliers):
    """Excluded contigs get a reason, and warnings are bulleted."""
    report = generate_extension_report(
        stats_dict,
        {},
        {'left_outliers': ['contig09'], 'right_outliers': []},
        {},
        ['chrM'],
        ['terminal screening window exceeds contig length'],
    )
    assert '## Excluded Contigs' in report
    assert '| chrM' in report
    assert 'user exclusion list' in report
    assert '## Warnings' in report
    assert 'terminal screening window exceeds contig length' in report


def test_anomalous_coverage_reported_separately_from_exclusions(stats_dict):
    """Flagged contigs are reported for review, not listed as excluded.

    Anomalous overhang coverage used to remove a contig from extension
    silently, and the report was passed empty outlier lists regardless, so the
    exclusion was never visible anywhere.
    """
    report = generate_extension_report(
        stats_dict,
        {},
        {'left_outliers': ['contig09'], 'right_outliers': []},
        {},
        ['chrM'],
        [],
    )

    assert '## Contigs With Anomalous Overhang Coverage' in report
    assert 'contig09' in report
    assert 'have **not** been excluded' in report

    # The flagged contig must not appear in the exclusion table.
    excluded_section = report.split('## Excluded Contigs')[1].split('##')[0]
    assert 'chrM' in excluded_section
    assert 'contig09' not in excluded_section


def test_anomaly_section_omitted_when_nothing_flagged(stats_dict, empty_outliers):
    """No anomaly section is emitted for an assembly with even coverage."""
    report = generate_extension_report(
        stats_dict,
        {},
        empty_outliers,
        {},
        [],
        [],
    )

    assert '## Contigs With Anomalous Overhang Coverage' not in report


def test_overhang_statistics_table(stats_dict, empty_outliers):
    """Overhang counts come from the contig stats, min/max render as integers."""
    overall = {
        'left': {'mean': 55.0, 'median': 55.0, 'std_dev': 7.07, 'min': 50, 'max': 60},
        'right': {'mean': 80.0, 'median': 80.0, 'std_dev': 0.0, 'min': 80, 'max': 80},
        'combined': {
            'mean': 63.3,
            'median': 60.0,
            'std_dev': 15.3,
            'min': 50,
            'max': 80,
        },
    }
    report = generate_extension_report(
        stats_dict,
        {},
        empty_outliers,
        overall,
        [],
        [],
    )
    assert '## Overhang Length Statistics' in report
    left_row = find_row(report, 'Left')
    assert left_row[1] == '2'  # two left overhangs in the fixture
    # Min/Max render as integers, not 50.00 / 60.00
    assert left_row[5] == '50'
    assert left_row[6] == '60'
    assert find_row(report, 'Combined')[1] == '3'


def test_per_contig_support_flags_extension_state(
    stats_dict, both_ends_extension, empty_outliers
):
    """The support table records which ends were actually extended."""
    report = generate_extension_report(
        stats_dict,
        both_ends_extension,
        empty_outliers,
        {},
        [],
        [],
    )
    row = find_row(report[report.index('## Per-Contig') :], 'contig01')
    assert row[2] == '2'  # left read count
    assert row[3] == '1'  # right read count
    assert row[4] == '60'  # longest left overhang
    assert row[5] == '80'  # longest right overhang
    assert row[6] == 'both'


def test_per_contig_support_marks_unextended(stats_dict, empty_outliers):
    """Contigs with support but no extension are marked 'no'."""
    report = generate_extension_report(
        stats_dict,
        {},
        empty_outliers,
        {},
        [],
        [],
    )
    row = find_row(report[report.index('## Per-Contig') :], 'contig01')
    assert row[6] == 'no'
