"""
Markdown statistics report for ``teloclip extend``.

The HTML report shows the evidence; this one is the record. It is plain
Markdown that also reads as aligned plain text in a terminal, so it can be
diffed between runs, pasted into a lab notebook, or grepped, and it states
every threshold that was in force so a result can be reproduced.

Table and number formatting primitives come from :mod:`teloclip.report.text`.
"""

from typing import Dict, List, Optional, Tuple

from ..core.analysis import ContigStats
from .text import fmt_delta, fmt_float, fmt_int, kv_table, md_table

# The HTML alignment view walks each read's CIGAR against the contig, so it
# needs real read bases rather than a pre-extracted anchor. These bound how much
# is retained and rendered; read slices are dropped as soon as a contig's panels
# are built, so nothing accumulates across the assembly.
#
# The rendered anchor window is at least HTML_MIN_WINDOW, or the longest anchor
# among the reads at that end if longer, capped at HTML_WINDOW_CAP and never
# more than the contig itself.
HTML_MIN_WINDOW = 500
# Generous enough to cover a typical long-read anchor whole; a pathological
# 50 kb anchor would otherwise render 50 kb of columns.
HTML_WINDOW_CAP = 2000
HTML_CONTIG_CONTEXT = HTML_WINDOW_CAP

# Read bases retained per overhang: enough for the widest window plus the clip
# and slack for insertions.
HTML_READ_CONTEXT = HTML_WINDOW_CAP + 1500

# Clipped bases rendered per read. Long telomeric overhangs run to kilobases,
# which no one reads base by base.
HTML_MAX_OVERHANG = 300


# Column order for the --overhang-log TSV.
OVERHANG_LOG_HEADER = (
    'contig',
    'contig_length',
    'end',
    'read',
    'aln_start',
    'aln_end',
    'gap_from_end',
    'clip_length',
    'overhang_length',
    'anchor_length',
)


def prune_read_slices(contig_stats: ContigStats, keep: int) -> None:
    """
    Drop retained read sequence from overhangs that will not be rendered.

    The HTML report shows at most ``keep`` reads per contig end, chosen by net
    gain. Every other overhang's slice is dead weight, and holding all of them
    across an assembly is what would make the report expensive on a real
    genome.

    Parameters
    ----------
    contig_stats : ContigStats
        Statistics for one contig. Modified in place.
    keep : int
        Reads per end whose sequence is retained.
    """
    for overhangs in (contig_stats.left_overhangs, contig_stats.right_overhangs):
        if len(overhangs) <= keep:
            continue
        ranked = sorted(
            overhangs, key=lambda oh: (oh.net_gain, oh.length), reverse=True
        )
        for oh in ranked[keep:]:
            oh.read_seq = ''
            oh.read_seq_offset = 0


def write_overhang_log_rows(handle, contig_stats: ContigStats) -> None:
    """
    Append one TSV row per accepted overhang for a contig.

    Written as the contig is processed rather than collected first, so the log
    costs no additional memory on large assemblies.

    Parameters
    ----------
    handle : file-like
        Open text handle positioned after the header row.
    contig_stats : ContigStats
        Overhang statistics for one contig.
    """
    for overhangs, end in (
        (contig_stats.left_overhangs, 'L'),
        (contig_stats.right_overhangs, 'R'),
    ):
        for oh in overhangs:
            # The gap is whatever part of the clip did not clear the terminus,
            # and is exactly the number of contig bases trimmed if used.
            gap = oh.length - oh.net_gain
            handle.write(
                '\t'.join(
                    str(field)
                    for field in (
                        contig_stats.contig_name,
                        contig_stats.contig_length,
                        end,
                        oh.read_name,
                        oh.alignment_pos,
                        oh.alignment_end,
                        gap,
                        oh.clip_length,
                        oh.net_gain,
                        oh.anchor_length,
                    )
                )
                + '\n'
            )


def _extension_lengths(ext_info: Dict) -> Tuple[int, int]:
    """
    Extract the net bases gained at each end from an extension record.

    Extending an end trims the contig bases the supporting read did not cover
    before grafting on its whole soft clip, so the bases the contig actually
    gains are the clip length minus the trim. Reporting the raw clip length
    instead overstates the extension and leaves the arithmetic in the report
    unable to reconcile with the final contig length.

    Parameters
    ----------
    ext_info : Dict
        Extension info dictionary as stored in ``extensions_applied``.

    Returns
    -------
    Tuple[int, int]
        Net bases gained at the left end and at the right end.
    """

    def net(side: str) -> int:
        """
        Net gain at one end, zero if that end was not extended.

        Parameters
        ----------
        side : str
            Either ``'left'`` or ``'right'``.

        Returns
        -------
        int
            Bases gained at that end.
        """
        if not ext_info.get(f'has_{side}_extension', False):
            return 0
        # apply_contig_extension records net_gain directly; fall back to the
        # clip-minus-trim arithmetic for dry runs and older records.
        if f'{side}_net_gain' in ext_info:
            return ext_info[f'{side}_net_gain']
        return ext_info.get(f'{side}_overhang_length', 0) - ext_info.get(
            f'{side}_trim_length', 0
        )

    return net('left'), net('right')


def _motif_gain_rows(
    extensions_applied: Dict[str, dict],
    terminal_motif_counts: Dict[str, Dict[str, Dict[str, int]]],
    post_motif_counts: Dict[str, Dict[str, Dict[str, int]]],
    screen_terminal_bases: int,
) -> Tuple[List[List[str]], int]:
    """
    Build the per-end motif analysis table rows.

    Parameters
    ----------
    extensions_applied : Dict[str, dict]
        Extension records keyed by contig name.
    terminal_motif_counts : Dict[str, Dict[str, Dict[str, int]]]
        Pre-extension motif counts, ``contig -> end -> motif -> count``, counted
        over the ``--screen-terminal-bases`` window.
    post_motif_counts : Dict[str, Dict[str, Dict[str, int]]]
        Post-extension motif counts over the same window plus the length of the
        extension at that end.
    screen_terminal_bases : int
        Size of the pre-extension screening window, in bases.

    Returns
    -------
    Tuple[List[List[str]], int]
        Table rows and the total motif gain summed across all contig ends.
    """
    rows: List[List[str]] = []
    total_gain = 0

    for contig_name in sorted(post_motif_counts):
        end_counts = post_motif_counts[contig_name]
        pre_counts = terminal_motif_counts.get(contig_name, {})
        left_added, right_added = _extension_lengths(
            extensions_applied.get(contig_name, {})
        )

        # The screening window is halved for contigs shorter than twice the
        # requested length, so report the window actually used rather than the
        # one asked for.
        effective_terminal = pre_counts.get('window', screen_terminal_bases)

        for end, added in (('left', left_added), ('right', right_added)):
            window = effective_terminal + added
            for motif_name in sorted(end_counts.get(end, {})):
                post = end_counts[end][motif_name]
                pre = pre_counts.get(end, {}).get(motif_name, 0)
                gain = post - pre
                total_gain += gain
                rows.append(
                    [
                        contig_name,
                        end,
                        motif_name,
                        fmt_int(window),
                        fmt_int(pre),
                        fmt_int(post),
                        fmt_delta(gain),
                    ]
                )

    return rows, total_gain


def generate_extension_report(
    stats_dict: Dict[str, ContigStats],
    extensions_applied: Dict[str, dict],
    outliers: Dict[str, List[str]],
    overall_stats: Dict[str, Dict[str, float]],
    excluded_contigs: List[str],
    warnings: List[str],
    motif_stats: Optional[Dict[str, Dict[str, int]]] = None,
    terminal_motif_counts: Optional[Dict[str, Dict[str, Dict[str, int]]]] = None,
    dry_run: bool = False,
    post_motif_counts: Optional[Dict[str, Dict[str, Dict[str, int]]]] = None,
    screen_terminal_bases: int = 0,
    total_contigs: Optional[int] = None,
) -> str:
    """
    Generate a Markdown statistics report for a contig extension run.

    The report opens with an at-a-glance summary and then presents per-contig
    detail as GitHub-flavoured Markdown tables, which also read as aligned
    plain text in a terminal.

    Parameters
    ----------
    stats_dict : Dict[str, ContigStats]
        Overhang statistics for every contig with overhang support.
    extensions_applied : Dict[str, dict]
        Extension records keyed by contig name.
    outliers : Dict[str, List[str]]
        Contigs flagged as having anomalous overhang coverage, under the keys
        ``left_outliers`` and ``right_outliers``. These are reported for review,
        not excluded.
    overall_stats : Dict[str, Dict[str, float]]
        Aggregate overhang statistics keyed by ``left``, ``right`` and
        ``combined``.
    excluded_contigs : List[str]
        Contigs excluded from extension by the user.
    warnings : List[str]
        Warning messages accumulated during the run.
    motif_stats : Optional[Dict[str, Dict[str, int]]], optional
        Whole-sequence motif counts per extended contig.
    terminal_motif_counts : Optional[Dict[str, Dict[str, Dict[str, int]]]], optional
        Pre-extension motif counts per contig end.
    dry_run : bool, optional
        Whether extensions were simulated rather than applied (default: False).
    post_motif_counts : Optional[Dict[str, Dict[str, Dict[str, int]]]], optional
        Post-extension motif counts per contig end, counted over the screening
        window plus the extension length at that end.
    screen_terminal_bases : int, optional
        Size of the terminal screening window in bases (default: 0).
    total_contigs : Optional[int], optional
        Number of contigs in the input assembly. Falls back to the number of
        contigs with overhang support when not supplied.

    Returns
    -------
    str
        The rendered Markdown report.
    """
    motif_stats = motif_stats or {}
    terminal_motif_counts = terminal_motif_counts or {}
    post_motif_counts = post_motif_counts or {}

    lines: List[str] = []

    def section(title: str, body: str) -> None:
        """
        Append a titled section to the report, skipping empty bodies.

        Parameters
        ----------
        title : str
            Section heading text, without the leading ``##``.
        body : str
            Rendered section body.

        Returns
        -------
        None
            The section is appended to the enclosing ``lines`` list.
        """
        if not body:
            return
        lines.append(f'## {title}')
        lines.append('')
        lines.append(body)
        lines.append('')

    # ------------------------------------------------------------------
    # Title
    # ------------------------------------------------------------------
    title = 'Teloclip Extend Report'
    if dry_run:
        title += ' (dry run)'
    lines.append(f'# {title}')
    lines.append('')

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    n_contigs = total_contigs if total_contigs is not None else len(stats_dict)
    ends_extended = 0
    total_gained = 0
    total_trimmed = 0

    for ext_info in extensions_applied.values():
        left_added, right_added = _extension_lengths(ext_info)
        ends_extended += int(left_added > 0) + int(right_added > 0)
        total_gained += left_added + right_added
        total_trimmed += ext_info.get('left_trim_length', 0) + ext_info.get(
            'right_trim_length', 0
        )

    motif_rows, total_motif_gain = _motif_gain_rows(
        extensions_applied,
        terminal_motif_counts,
        post_motif_counts,
        screen_terminal_bases,
    )

    motif_names = sorted(
        {
            motif
            for counts in post_motif_counts.values()
            for end in counts.values()
            for motif in end
        }
    )

    summary_pairs = [
        ('Mode', 'dry run (no sequences written)' if dry_run else 'extension applied'),
        ('Contigs in assembly', fmt_int(n_contigs)),
        ('Contigs with overhang support', fmt_int(len(stats_dict))),
        ('Contigs extended', fmt_int(len(extensions_applied))),
        (
            'Contig ends extended',
            f'{fmt_int(ends_extended)} of {fmt_int(n_contigs * 2)}',
        ),
        ('Net bases gained', fmt_int(total_gained)),
        ('Bases trimmed back', fmt_int(total_trimmed)),
        ('Raw overhang bases grafted', fmt_int(total_gained + total_trimmed)),
    ]
    if motif_names:
        summary_pairs.append(('Motifs counted', ', '.join(motif_names)))
        summary_pairs.append(('Total motif gain', fmt_delta(total_motif_gain)))
    if screen_terminal_bases > 0:
        summary_pairs.append(
            ('Terminal screening window', f'{fmt_int(screen_terminal_bases)} bp')
        )
    summary_pairs.append(('Contigs excluded', fmt_int(len(excluded_contigs))))
    summary_pairs.append(('Warnings', fmt_int(len(warnings))))

    section('Summary', kv_table(summary_pairs))

    # ------------------------------------------------------------------
    # Extensions applied
    # ------------------------------------------------------------------
    ext_rows: List[List[str]] = []
    for contig_name in sorted(extensions_applied):
        ext_info = extensions_applied[contig_name]
        left_added, right_added = _extension_lengths(ext_info)
        ext_rows.append(
            [
                contig_name,
                fmt_int(ext_info.get('original_length', 0)),
                fmt_int(ext_info.get('final_length', 0)),
                fmt_delta(left_added),
                fmt_delta(right_added),
                fmt_delta(left_added + right_added),
                ext_info.get('left_read_name', '-') if left_added else '-',
                ext_info.get('right_read_name', '-') if right_added else '-',
                fmt_int(ext_info.get('left_trim_length', 0)),
                fmt_int(ext_info.get('right_trim_length', 0)),
            ]
        )

    ext_note = (
        'The `+bp` columns are net: the overhang grafted on at that end, less '
        'the contig bases trimmed to make room for it. So '
        '`Original bp + Total +bp = Final bp` for every row.\n\n'
        if ext_rows
        else ''
    )

    section(
        'Extensions That Would Be Applied' if dry_run else 'Extensions Applied',
        ext_note
        + md_table(
            [
                'Contig',
                'Original bp',
                'Final bp',
                'Left +bp',
                'Right +bp',
                'Total +bp',
                'Left read',
                'Right read',
                'Left trim',
                'Right trim',
            ],
            ext_rows,
            align=['l', 'r', 'r', 'r', 'r', 'r', 'l', 'l', 'r', 'r'],
        ),
    )

    # ------------------------------------------------------------------
    # Motif analysis
    # ------------------------------------------------------------------
    if motif_rows:
        lines.append('## Telomere Motif Analysis')
        lines.append('')
        lines.append(
            f'Counts are per contig end. `Pre` counts the {fmt_int(screen_terminal_bases)} bp '
            'terminal screening window of the original contig; `Post` counts that same window '
            'plus the bases added at that end, in the extended contig.'
        )
        lines.append('')
        lines.append(
            md_table(
                ['Contig', 'End', 'Motif', 'Window bp', 'Pre', 'Post', 'Gain'],
                motif_rows,
                align=['l', 'l', 'l', 'r', 'r', 'r', 'r'],
            )
        )
        lines.append('')
        lines.append(f'**Total motif gain: {fmt_delta(total_motif_gain)}**')
        lines.append('')

    # ------------------------------------------------------------------
    # Overhang statistics
    # ------------------------------------------------------------------
    n_left = sum(cs.left_count for cs in stats_dict.values())
    n_right = sum(cs.right_count for cs in stats_dict.values())
    set_counts = {'left': n_left, 'right': n_right, 'combined': n_left + n_right}

    stat_rows: List[List[str]] = []
    for label in ('left', 'right', 'combined'):
        stats = overall_stats.get(label)
        if not stats:
            continue
        stat_rows.append(
            [
                label.title(),
                fmt_int(set_counts[label]),
                fmt_float(stats.get('mean', 0.0)),
                fmt_float(stats.get('median', 0.0)),
                fmt_float(stats.get('std_dev', 0.0)),
                fmt_int(int(stats.get('min', 0))),
                fmt_int(int(stats.get('max', 0))),
            ]
        )

    section(
        'Overhang Length Statistics',
        md_table(
            ['Set', 'N', 'Mean', 'Median', 'SD', 'Min', 'Max'],
            stat_rows,
            align=['l', 'r', 'r', 'r', 'r', 'r', 'r'],
        ),
    )

    # ------------------------------------------------------------------
    # Per-contig overhang support
    # ------------------------------------------------------------------
    support_rows: List[List[str]] = []
    for contig_name in sorted(stats_dict):
        contig_stats = stats_dict[contig_name]
        longest_left = max((oh.length for oh in contig_stats.left_overhangs), default=0)
        longest_right = max(
            (oh.length for oh in contig_stats.right_overhangs), default=0
        )
        ext_info = extensions_applied.get(contig_name, {})
        left_added, right_added = _extension_lengths(ext_info)
        if left_added and right_added:
            extended = 'both'
        elif left_added:
            extended = 'left'
        elif right_added:
            extended = 'right'
        else:
            extended = 'no'
        support_rows.append(
            [
                contig_name,
                fmt_int(contig_stats.contig_length),
                fmt_int(contig_stats.left_count),
                fmt_int(contig_stats.right_count),
                fmt_int(longest_left),
                fmt_int(longest_right),
                extended,
            ]
        )

    if support_rows:
        lines.append('## Per-Contig Overhang Support')
        lines.append('')
        lines.append(
            'Read counts are clipped reads anchored at that contig end which passed '
            'filtering. `Longest` is the longest overhang seen at that end, which is the '
            'sequence used for extension unless it failed a homopolymer or length check.'
        )
        lines.append('')
        lines.append(
            md_table(
                [
                    'Contig',
                    'Length',
                    'Left reads',
                    'Right reads',
                    'Longest left',
                    'Longest right',
                    'Extended',
                ],
                support_rows,
                align=['l', 'r', 'r', 'r', 'r', 'r', 'l'],
            )
        )
        lines.append('')

    # ------------------------------------------------------------------
    # Excluded contigs
    # ------------------------------------------------------------------
    section(
        'Excluded Contigs',
        md_table(
            ['Contig', 'Reason'],
            [[contig, 'user exclusion list'] for contig in excluded_contigs],
            align=['l', 'l'],
        ),
    )

    # ------------------------------------------------------------------
    # Anomalous overhang coverage
    # ------------------------------------------------------------------
    anomaly_rows: List[List[str]] = []
    for end in ('left', 'right'):
        for contig in outliers.get(f'{end}_outliers', []):
            counts = stats_dict.get(contig)
            observed = getattr(counts, f'{end}_count', 0) if counts is not None else 0
            anomaly_rows.append([contig, end, fmt_int(observed)])

    if anomaly_rows:
        section(
            'Contigs With Anomalous Overhang Coverage',
            'These contig ends carry far more clipped reads than the rest of '
            'the assembly. That often marks a collapsed repeat, an rDNA array '
            'or an organellar contig attracting reads from elsewhere, in which '
            'case extending from a single read is not meaningful.\n\n'
            'They have **not** been excluded. Review them and, if you agree, '
            're-run with `--exclude-contigs`.\n\n'
            + md_table(
                ['Contig', 'End', 'Overhang reads'],
                anomaly_rows,
                align=['l', 'l', 'r'],
            ),
        )

    # ------------------------------------------------------------------
    # Warnings
    # ------------------------------------------------------------------
    if warnings:
        lines.append('## Warnings')
        lines.append('')
        for warning in warnings:
            lines.append(f'- {warning}')
        lines.append('')

    if extensions_applied and not dry_run:
        lines.append(
            'Note: extended contigs should be polished (e.g. Medaka for ONT data, '
            'Pypolca for Illumina data) before downstream analysis.'
        )
        lines.append('')

    return '\n'.join(lines).rstrip() + '\n'
