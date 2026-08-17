"""Unit tests for the teloclip.report.text formatting helpers."""

import pytest

from teloclip.report.text import fmt_delta, fmt_float, fmt_int, kv_table, md_table


def test_fmt_int_adds_thousands_separators():
    """Integers are rendered with comma separators."""
    assert fmt_int(0) == '0'
    assert fmt_int(999) == '999'
    assert fmt_int(45120) == '45,120'
    assert fmt_int(6031340) == '6,031,340'


def test_fmt_delta_signs_and_zero():
    """Deltas are signed, and zero renders as a dash."""
    assert fmt_delta(180) == '+180'
    assert fmt_delta(-12) == '-12'
    assert fmt_delta(0) == '-'
    assert fmt_delta(1234) == '+1,234'


def test_fmt_float_decimal_places():
    """Floats honour the requested number of decimal places."""
    assert fmt_float(58.888) == '58.9'
    assert fmt_float(58.888, dp=2) == '58.89'
    assert fmt_float(1234.5) == '1,234.5'


def test_md_table_empty_rows_returns_empty_string():
    """An empty body yields an empty string so callers can skip the section."""
    assert md_table(['A', 'B'], []) == ''


def test_md_table_separator_encodes_alignment():
    """The separator row uses the GFM markers for each alignment code."""
    table = md_table(
        ['Left', 'Right', 'Centre'],
        [['a', 'b', 'c']],
        align=['l', 'r', 'c'],
    )
    separator = table.split('\n')[1]
    cells = [cell.strip() for cell in separator.strip('|').split('|')]
    assert cells == ['---', '---:', ':---:']


def test_md_table_pads_columns_to_common_width():
    """Every rendered line is the same width, so raw text stays aligned."""
    table = md_table(
        ['Contig', 'Length'],
        [['contig01', '560'], ['a_much_longer_contig', '1,610']],
        align=['l', 'r'],
    )
    lines = table.split('\n')
    assert len({len(line) for line in lines}) == 1
    # Right-aligned values are flushed to the right of their column.
    assert lines[2].endswith('|    560 |')


def test_md_table_row_width_mismatch_raises():
    """A row with the wrong number of cells is an error, not silent truncation."""
    with pytest.raises(ValueError, match='2 cells'):
        md_table(['A', 'B', 'C'], [['x', 'y']])


def test_md_table_bad_align_code_raises():
    """Unknown alignment codes are rejected."""
    with pytest.raises(ValueError, match='alignment code'):
        md_table(['A'], [['x']], align=['x'])


def test_md_table_align_length_mismatch_raises():
    """The align sequence must match the number of columns."""
    with pytest.raises(ValueError, match='2 entries'):
        md_table(['A', 'B', 'C'], [['x', 'y', 'z']], align=['l', 'r'])


def test_kv_table_shape():
    """kv_table renders a two-column Metric/Value table."""
    table = kv_table([('Contigs extended', '4'), ('Total bases added', '551')])
    lines = table.split('\n')
    assert lines[0].startswith('| Metric')
    assert 'Value' in lines[0]
    assert len(lines) == 4  # header, separator, two rows
    assert 'Contigs extended' in lines[2]


def test_kv_table_custom_key_header():
    """The label column header can be overridden."""
    table = kv_table([('a', 'b')], key_header='Setting')
    assert table.split('\n')[0].startswith('| Setting')
