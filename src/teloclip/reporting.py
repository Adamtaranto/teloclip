"""
Plain-text and Markdown formatting helpers for teloclip reports.

This module provides small, dependency-free helpers used to build the
statistics reports emitted by the teloclip sub-commands. Tables are written as
GitHub-flavoured Markdown, padded so that they also read as aligned plain text
in a terminal.
"""

from typing import List, Optional, Sequence, Tuple

# Separator cell contents for each supported alignment code.
_ALIGN_MARKERS = {
    'l': '---',
    'r': '---:',
    'c': ':---:',
}


def fmt_int(value: int) -> str:
    """
    Format an integer with thousands separators.

    Parameters
    ----------
    value : int
        Value to format.

    Returns
    -------
    str
        The value rendered with comma thousands separators, e.g. ``'45,120'``.
    """
    return f'{value:,}'


def fmt_delta(value: int) -> str:
    """
    Format a signed change in length.

    Parameters
    ----------
    value : int
        Change in length, in base pairs.

    Returns
    -------
    str
        A signed, comma-separated string such as ``'+180'`` or ``'-12'``.
        Zero is rendered as ``'-'`` so that unchanged ends are visually quiet.
    """
    if value == 0:
        return '-'
    return f'{value:+,}'


def fmt_float(value: float, dp: int = 1) -> str:
    """
    Format a float to a fixed number of decimal places.

    Parameters
    ----------
    value : float
        Value to format.
    dp : int, optional
        Number of decimal places to show (default: 1).

    Returns
    -------
    str
        The value rendered with ``dp`` decimal places.
    """
    return f'{value:,.{dp}f}'


def md_table(
    headers: Sequence[str],
    rows: Sequence[Sequence[str]],
    align: Optional[Sequence[str]] = None,
) -> str:
    """
    Render a GitHub-flavoured Markdown table.

    Cells are padded to a common width per column so that the raw text is also
    aligned when viewed in a terminal.

    Parameters
    ----------
    headers : Sequence[str]
        Column header labels.
    rows : Sequence[Sequence[str]]
        Table body. Each row must have the same number of cells as ``headers``.
    align : Optional[Sequence[str]], optional
        Per-column alignment codes, one of ``'l'``, ``'r'`` or ``'c'``. Defaults
        to left alignment for every column.

    Returns
    -------
    str
        The rendered table, or an empty string when ``rows`` is empty so that
        callers can omit the surrounding section.

    Raises
    ------
    ValueError
        If a row or the ``align`` sequence does not match the header width, or
        if an unknown alignment code is given.
    """
    if not rows:
        return ''

    ncols = len(headers)

    if align is None:
        align = ['l'] * ncols
    elif len(align) != ncols:
        raise ValueError(
            f'align has {len(align)} entries but there are {ncols} columns'
        )

    for code in align:
        if code not in _ALIGN_MARKERS:
            raise ValueError(f'Unknown alignment code: {code!r}')

    cells: List[List[str]] = [[str(h) for h in headers]]
    for row in rows:
        if len(row) != ncols:
            raise ValueError(
                f'Row has {len(row)} cells but there are {ncols} columns: {row!r}'
            )
        cells.append([str(c) for c in row])

    widths = [max(len(row[i]) for row in cells) for i in range(ncols)]
    # A right-aligned separator needs at least 4 characters ('---:').
    widths = [max(w, len(_ALIGN_MARKERS[align[i]])) for i, w in enumerate(widths)]

    def _render(row: Sequence[str]) -> str:
        """
        Render one row with each cell padded to its column width.

        Parameters
        ----------
        row : Sequence[str]
            Cell values for a single row.

        Returns
        -------
        str
            The pipe-delimited, padded row.
        """
        padded = []
        for i, cell in enumerate(row):
            if align[i] == 'r':
                padded.append(cell.rjust(widths[i]))
            elif align[i] == 'c':
                padded.append(cell.center(widths[i]))
            else:
                padded.append(cell.ljust(widths[i]))
        return '| ' + ' | '.join(padded) + ' |'

    lines = [_render(cells[0])]
    lines.append(
        '| '
        + ' | '.join(
            _ALIGN_MARKERS[align[i]].rjust(widths[i])
            if align[i] == 'r'
            else _ALIGN_MARKERS[align[i]].ljust(widths[i])
            for i in range(ncols)
        )
        + ' |'
    )
    lines.extend(_render(row) for row in cells[1:])

    return '\n'.join(lines)


def kv_table(pairs: Sequence[Tuple[str, str]], key_header: str = 'Metric') -> str:
    """
    Render a two-column metric/value table.

    Parameters
    ----------
    pairs : Sequence[Tuple[str, str]]
        Ordered ``(label, value)`` pairs.
    key_header : str, optional
        Header for the label column (default: ``'Metric'``).

    Returns
    -------
    str
        The rendered table, or an empty string when ``pairs`` is empty.
    """
    return md_table(
        [key_header, 'Value'],
        [[label, value] for label, value in pairs],
        align=['l', 'r'],
    )


def histogram(
    counts: Sequence[int],
    bins: int = 10,
    width: int = 24,
    label: str = 'count',
) -> str:
    """
    Render a horizontal ASCII histogram of a distribution.

    Intended for terminal logs, where a per-contig distribution is far easier to
    read as a shape than as a list of numbers. Returns plain text with no colour
    codes so it survives redirection to a log file.

    Parameters
    ----------
    counts : Sequence[int]
        Observed values to bin.
    bins : int, optional
        Number of bins (default: 10). Reduced automatically when the observed
        range spans fewer distinct integers.
    width : int, optional
        Maximum bar width in characters (default: 24, chosen so the result
        still fits an 80-column terminal alongside a rich log prefix).
    label : str, optional
        Noun describing what is being counted, used in the header
        (default: ``'count'``).

    Returns
    -------
    str
        The rendered histogram, or an empty string when ``counts`` is empty.
    """
    values = list(counts)
    if not values:
        return ''

    low, high = min(values), max(values)

    # A distribution with no spread is a single bar; binning it would produce
    # a row of empty ranges.
    if low == high:
        return f'{label} {low}: {"#" * min(width, len(values))} ({len(values)})'

    # Never use more bins than there are distinct integer values to put in them.
    bins = max(1, min(bins, high - low + 1))
    span = (high - low + 1) / bins

    binned = [0] * bins
    for value in values:
        index = min(bins - 1, int((value - low) / span))
        binned[index] += 1

    peak = max(binned)
    edges = [(low + int(i * span), low + int((i + 1) * span) - 1) for i in range(bins)]
    edge_width = max(len(f'{lo}-{hi}') for lo, hi in edges)
    count_width = len(str(peak))

    lines = []
    for (lo, hi), n in zip(edges, binned):
        bar = '#' * int(round(width * n / peak)) if peak else ''
        lines.append(
            f'  {f"{lo}-{hi}":>{edge_width}} | {bar:<{width}} {n:>{count_width}}'
        )

    return '\n'.join(lines)
