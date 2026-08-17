"""
Extend sub-command implementation.

This module implements the 'teloclip extend' command for automatically extending
draft contigs using overhang analysis from soft-clipped alignments.

This version is optimized for large genomes using streaming I/O and indexed access
to avoid loading entire genomes into memory.
"""

from datetime import datetime
import logging
from pathlib import Path
import re
import shlex
import sys
from typing import Dict, List, Optional, Tuple

import click
import pyfaidx
import pysam

from ..core.analysis import (
    ContigStats,
    calculate_overhang_statistics,
    flag_anomalous_overhang_coverage,
)
from ..core.motifs import make_fuzzy_motif_regex, make_motif_regex
from ..core.seqops import read_fai, revComp
from ..core.streaming_analysis import (
    process_single_contig_extension,
    stream_contigs_for_extension,
)
from ..io.streaming import (
    BufferedContigWriter,
    StreamingGenomeProcessor,
    validate_indexed_files,
)
from ..report.text import fmt_delta, fmt_float, fmt_int, histogram, kv_table, md_table


def validate_output_directories(output_fasta: Path, stats_report: Path) -> None:
    """
    Validate that output directories exist and are writable.

    Parameters
    ----------
    output_fasta : Path
        Path where extended FASTA will be written.
    stats_report : Path
        Path where statistics report will be written.

    Returns
    -------
    None
        Function validates and creates output directories as needed.

    Raises
    ------
    click.ClickException
        If output directories cannot be created or are not writable.
    """
    for output_path in [output_fasta, stats_report]:
        if output_path:
            output_dir = output_path.parent
            if not output_dir.exists():
                try:
                    output_dir.mkdir(parents=True, exist_ok=True)
                except OSError as e:
                    raise click.ClickException(
                        f'Cannot create output directory {output_dir}: {e}'
                    ) from e


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


def _prune_read_slices(contig_stats: ContigStats, keep: int) -> None:
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


def _write_overhang_log_rows(handle, contig_stats: ContigStats) -> None:
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


def count_terminal_motifs(
    fasta_file: Path,
    contig_dict: Dict[str, int],
    motif_patterns: Dict[str, re.Pattern],
    terminal_length: int,
) -> Dict[str, Dict[str, Dict[str, int]]]:
    """
    Count motifs in terminal regions of contigs.

    Parameters
    ----------
    fasta_file : Path
        Path to indexed FASTA file.
    contig_dict : Dict[str, int]
        Dictionary mapping contig names to their lengths.
    motif_patterns : Dict[str, re.Pattern]
        Dictionary of compiled motif regex patterns.
    terminal_length : int
        Number of bases to extract from each terminal end.

    Returns
    -------
    Dict[str, Dict[str, Dict[str, int]]]
        Nested dictionary: contig_name -> {left: {motif: count}, right: {motif: count}}.

    Raises
    ------
    click.ClickException
        If FASTA file cannot be opened or contigs cannot be accessed.
    """
    terminal_counts = {}
    warnings_issued = []

    if terminal_length <= 0 or not motif_patterns:
        return terminal_counts

    try:
        # Open FASTA file with pyfaidx
        fasta = pyfaidx.Fasta(str(fasta_file))

        for contig_name, contig_length in contig_dict.items():
            # Check if terminal length is too large relative to contig
            if terminal_length > contig_length * 0.5:
                warning_msg = (
                    f'Terminal screening length ({terminal_length}) is > 50% '
                    f'of contig {contig_name} length ({contig_length})'
                )
                if warning_msg not in warnings_issued:
                    logging.warning(warning_msg)
                    warnings_issued.append(warning_msg)

            # Adjust terminal length for short contigs
            actual_length = min(terminal_length, contig_length // 2)
            if actual_length <= 0:
                continue

            try:
                # Get terminal sequences
                left_seq = str(fasta[contig_name][:actual_length]).upper()
                right_seq = str(fasta[contig_name][-actual_length:]).upper()

                # Initialize counts for this contig. 'window' records the
                # window actually used, which is halved for short contigs and
                # so may be smaller than the requested terminal_length.
                terminal_counts[contig_name] = {
                    'left': {},
                    'right': {},
                    'window': actual_length,
                }

                # Count motifs in each terminal region
                for motif_name, pattern in motif_patterns.items():
                    left_matches = len(pattern.findall(left_seq))
                    right_matches = len(pattern.findall(right_seq))

                    terminal_counts[contig_name]['left'][motif_name] = left_matches
                    terminal_counts[contig_name]['right'][motif_name] = right_matches

            except (KeyError, IndexError) as e:
                logging.warning(
                    f'Could not access contig {contig_name} in FASTA file: {e}'
                )
                continue

        fasta.close()

    except Exception as e:
        raise click.ClickException(
            f'Error reading FASTA file for terminal motif screening: {e}'
        ) from e

    return terminal_counts


def read_excluded_contigs_file(exclude_file: Path) -> List[str]:
    """
    Read contig names from a file, handling different line endings.

    Parameters
    ----------
    exclude_file : Path
        Path to file containing contig names (one per line).

    Returns
    -------
    List[str]
        List of contig names from the file.

    Raises
    ------
    click.ClickException
        If file cannot be read or is empty.
    """
    try:
        # Read file with universal newlines to handle different line endings
        with open(exclude_file, 'r', newline=None) as f:
            raw_lines = f.readlines()

        # Process lines: strip whitespace and filter out empty lines
        contig_names = []
        for line_num, line in enumerate(raw_lines, 1):
            # Strip all whitespace (including \r, \n, spaces, tabs)
            cleaned_line = line.strip()

            # Skip empty lines
            if not cleaned_line:
                logging.debug(f'Skipping empty line {line_num} in {exclude_file}')
                continue

            contig_names.append(cleaned_line)

        if not contig_names:
            raise click.ClickException(
                f'No valid contig names found in exclusion file: {exclude_file}'
            )

        logging.info(f'Read {len(contig_names)} contig names from {exclude_file}')
        return contig_names

    except FileNotFoundError:
        raise click.ClickException(
            f'Exclusion file not found: {exclude_file}'
        ) from None
    except IOError as e:
        raise click.ClickException(
            f'Error reading exclusion file {exclude_file}: {e}'
        ) from e


def combine_excluded_contigs(
    exclude_contigs_str: str,
    exclude_contigs_file: Path,
    contig_dict: Dict[str, int],
) -> set:
    """
    Combine excluded contig names from string and file sources.

    Parameters
    ----------
    exclude_contigs_str : str
        Comma-delimited string of contig names to exclude.
    exclude_contigs_file : Path
        Path to file containing contig names to exclude.
    contig_dict : Dict[str, int]
        Dictionary mapping contig names to their lengths.

    Returns
    -------
    set
        Combined set of contig names to exclude.

    Raises
    ------
    click.ClickException
        If any listed name is absent from the reference index.
    """
    all_excluded_names = []

    # Collect from string source
    if exclude_contigs_str:
        string_contigs = [
            contig.strip()
            for contig in exclude_contigs_str.split(',')
            if contig.strip()
        ]
        if string_contigs:
            all_excluded_names.extend(string_contigs)
            logging.info(f'Found {len(string_contigs)} contigs from --exclude-contigs')

    # Collect from file source
    if exclude_contigs_file:
        file_contigs = read_excluded_contigs_file(exclude_contigs_file)
        if file_contigs:
            all_excluded_names.extend(file_contigs)
            logging.info(
                f'Found {len(file_contigs)} contigs from --exclude-contigs-file'
            )

    # Check if both sources provided and warn
    if exclude_contigs_str and exclude_contigs_file:
        logging.warning(
            'Both --exclude-contigs and --exclude-contigs-file provided. '
            'Combining contig names from both sources.'
        )

    if not all_excluded_names:
        return set()

    # Create unique set and validate against contig dictionary
    unique_names = set(all_excluded_names)
    duplicates_removed = len(all_excluded_names) - len(unique_names)
    if duplicates_removed > 0:
        logging.info(f'Removed {duplicates_removed} duplicate contig names')

    # Validate every name against the reference index. A misspelled name
    # excludes nothing while looking like it worked, so this is an error rather
    # than a warning: the run would otherwise extend a contig the user believed
    # they had held back.
    unknown = sorted(name for name in unique_names if name not in contig_dict)
    if unknown:
        raise click.ClickException(
            'These contigs were listed for exclusion but are not in the '
            f'reference index: {", ".join(unknown)}. '
            'Check the spelling against the .fai file.'
        )

    excluded_set = set(unique_names)
    for contig_name in sorted(excluded_set):
        logging.info(f'Contig "{contig_name}" will be excluded from extension')
    logging.info(f'Total unique contigs excluded: {len(excluded_set)}')

    return excluded_set


def get_motif_regex(motif_str: str, fuzzy: bool = False) -> Dict[str, re.Pattern]:
    """
    Generate motif regex patterns from a comma-delimited string.

    Parameters
    ----------
    motif_str : str
        Comma-delimited motif sequences.
    fuzzy : bool, optional
        Whether to use fuzzy matching allowing ±1 character variation. Default is False.

    Returns
    -------
    Dict[str, re.Pattern]
        Dictionary mapping motif sequences to compiled regex patterns.
    """
    # Initialize motif patterns dictionary
    motif_patterns = {}
    # Parse comma-delimited motifs
    raw_motifs = [
        motif.strip().upper() for motif in motif_str.split(',') if motif.strip()
    ]
    logging.info(f'Processing motif list: {", ".join(raw_motifs)}')

    # Validate motifs - must contain only A, T, G, C
    valid_bases = {'A', 'T', 'G', 'C'}
    validated_motifs = []

    for motif in raw_motifs:
        if not motif:  # Skip empty motifs
            continue
        if not all(base in valid_bases for base in motif):
            invalid_bases = set(motif) - valid_bases
            logging.warning(
                f'Skipping invalid motif "{motif}": contains invalid bases {invalid_bases}'
            )
            continue
        validated_motifs.append(motif)

    if not validated_motifs:
        logging.warning('No valid motifs found after validation')
    else:
        logging.info(
            f'Validated {len(validated_motifs)} motifs: {", ".join(validated_motifs)}'
        )

        # Add reverse complements and create unique set
        all_motifs = set()
        for motif in validated_motifs:
            all_motifs.add(motif)
            rev_comp = revComp(motif)
            all_motifs.add(rev_comp)
            logging.debug(f'Motif: {motif} -> Reverse complement: {rev_comp}')

        # Convert to sorted list for consistent ordering
        unique_motifs = sorted(all_motifs)
        logging.info(
            f'Final motif set (including reverse complements): {", ".join(unique_motifs)}'
        )

        # Create regex patterns for each unique motif
        for motif in unique_motifs:
            if fuzzy:
                pattern_str = make_fuzzy_motif_regex(motif)
                pattern_name = f'{motif}_fuzzy'
                logging.debug(f'Created fuzzy pattern for {motif}: {pattern_str}')
            else:
                pattern_str = make_motif_regex(motif)
                pattern_name = motif
                logging.debug(f'Created exact pattern for {motif}: {pattern_str}')

            # Compile the pattern for use with re.findall
            motif_patterns[pattern_name] = re.compile(pattern_str)

        logging.info(f'Created {len(motif_patterns)} motif patterns for analysis')
        if fuzzy:
            logging.info(
                'Using fuzzy matching (±1 character variation) for motif counting'
            )

    return motif_patterns


@click.command(
    help='Extend contigs using overhang analysis from soft-clipped alignments.'
)
@click.argument('bam_file', type=click.Path(exists=True, path_type=Path))
@click.argument('reference_fasta', type=click.Path(exists=True, path_type=Path))
@click.option(
    '--output-fasta', type=click.Path(path_type=Path), help='Extended FASTA output file'
)
@click.option(
    '--stats-report',
    type=click.Path(path_type=Path),
    help='Statistics report output file',
)
@click.option(
    '--exclude-outliers',
    is_flag=True,
    help='DEPRECATED and ignored. Contigs with anomalous overhang coverage are '
    'now reported for review rather than silently dropped; exclude them with '
    '--exclude-contigs if you agree with the assessment.',
)
@click.option(
    '--outlier-threshold',
    type=float,
    default=3.5,
    help='Modified z-score above which a contig end is reported as having '
    'anomalous overhang coverage (default: 3.5)',
)
@click.option(
    '--min-overhangs',
    type=int,
    default=1,
    help='Minimum supporting overhangs required (default: 1)',
)
@click.option(
    '--max-homopolymer',
    type=int,
    default=500,
    help='Maximum homopolymer run length allowed (default: 500)',
)
@click.option(
    '--min-extension',
    type=int,
    default=1,
    help='Minimum novel bases an overhang must contribute to be used (default: 1)',
)
@click.option(
    '--min-clip',
    type=int,
    default=1,
    help='Require clip to extend past the contig end by at least N bases (default: 1)',
)
@click.option(
    '--max-break',
    type=int,
    default=50,
    help='Maximum gap allowed between alignment and contig end (default: 50)',
)
@click.option(
    '--min-anchor',
    type=int,
    default=100,
    help='Minimum anchor length required for alignment (default: 100)',
)
@click.option(
    '--dry-run', is_flag=True, help='Report extensions without modifying sequences'
)
@click.option(
    '--count-motifs',
    type=str,
    help='Comma-delimited motif sequences to count in overhang regions (e.g., "TTAGGG,CCCTAA")',
)
@click.option(
    '--fuzzy-count',
    is_flag=True,
    help='Use fuzzy motif matching allowing ±1 character variation when counting motifs',
)
@click.option(
    '--prefix',
    type=str,
    default='teloclip_extended',
    help='Prefix for default output filenames (default: teloclip_extended)',
)
@click.option(
    '--screen-terminal-bases',
    type=int,
    default=0,
    help='Number of terminal bases to screen for motifs in original contigs (default: 0, disabled)',
)
@click.option(
    '--exclude-contigs',
    type=str,
    help='Comma-delimited list of contig names to exclude from extension (e.g., "chrM,chrC,scaffold_123")',
)
@click.option(
    '--exclude-contigs-file',
    type=click.Path(exists=True, path_type=Path),
    help='Text file containing contig names to exclude (one per line)',
)
@click.option(
    '--log-level',
    default='INFO',
    type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR'], case_sensitive=False),
    help='Logging level (default: INFO).',
)
@click.option(
    '--logfile',
    type=click.Path(path_type=Path),
    help='Also write log messages to this file (parent directories are created).',
)
@click.option(
    '--html-report',
    type=click.Path(path_type=Path),
    help='Write a self-contained HTML report showing every overhang read '
    'aligned against the contig terminus it supports, plus overhang depth '
    'across the assembly.',
)
@click.option(
    '--html-max-reads',
    type=int,
    default=25,
    help='Maximum overhang reads rendered per contig end in the HTML report '
    '(default: 25). Reads contributing the most sequence are shown first.',
)
@click.option(
    '--overhang-log',
    type=click.Path(path_type=Path),
    help='Write a TSV describing every accepted overhang read: contig, end, '
    'gap from the contig terminus, clip length and overhang length.',
)
@click.pass_context
def extend(
    ctx,
    bam_file,
    reference_fasta,
    output_fasta,
    stats_report,
    exclude_outliers,
    outlier_threshold,
    min_overhangs,
    max_homopolymer,
    min_extension,
    min_clip,
    max_break,
    min_anchor,
    dry_run,
    count_motifs,
    fuzzy_count,
    prefix,
    screen_terminal_bases,
    exclude_contigs,
    exclude_contigs_file,
    log_level,
    logfile,
    html_report,
    html_max_reads,
    overhang_log,
):
    """
    Extend contigs based on alignment overhangs.

    This command analyzes soft-clipped reads aligned to contig ends to extend
    sequences where there is sufficient evidence of consensus sequence beyond
    the current contig boundaries.

    Parameters
    ----------
    ctx : click.Context
        Click context object.
    bam_file : str
        Path to sorted and indexed BAM file.
    reference_fasta : str
        Path to reference FASTA file (must be indexed).
    output_fasta : str
        Path for output extended FASTA file.
    stats_report : str
        Path for output statistics report.
    exclude_outliers : bool
        Deprecated and ignored. Retained so existing command lines keep working;
        emits a warning pointing at --exclude-contigs.
    outlier_threshold : float
        Modified z-score above which a contig end is reported as having
        anomalous overhang coverage.
    min_overhangs : int
        Minimum number of overhangs required for extension.
    max_homopolymer : int
        Maximum homopolymer length allowed in extensions.
    min_extension : int
        Minimum number of novel bases an overhang must contribute to be used.
    min_clip : int
        Minimum number of clipped bases required past the contig terminus.
    max_break : int
        Maximum tolerated gap between a contig terminus and the alignment.
    min_anchor : int
        Minimum number of anchoring (M/=/X) bases required.
    dry_run : bool
        Perform analysis without writing output files.
    count_motifs : bool
        Count telomeric motifs in extensions.
    fuzzy_count : int
        Number of mismatches allowed in motif counting.
    prefix : str
        Prefix for default output filenames.
    screen_terminal_bases : int
        Number of terminal bases to screen for motifs in original contigs.
    exclude_contigs : str
        Comma-delimited list of contig names to exclude from extension.
    exclude_contigs_file : Path
        Path to file containing contig names to exclude (one per line).
    log_level : str
        Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
    logfile : Path or None
        Optional path to also write log messages to.
    html_report : Path or None
        Optional path for a self-contained HTML report.
    html_max_reads : int
        Maximum overhang reads rendered per contig end in the HTML report.
    overhang_log : Path or None
        Optional path for a per-overhang TSV describing every accepted read.
    """
    from ..logs import init_logging

    # Initialize logging for this command
    init_logging(log_level, logfile)

    ctx.ensure_object(dict)

    # A clip that does not reach past the contig terminus contributes no novel
    # sequence, and applying it would trim more bases than it adds. Clamping
    # here is what makes "every accepted overhang lengthens the contig" hold.
    if min_clip < 1:
        logging.warning(
            f'--min-clip must be at least 1 (got {min_clip}); using 1. '
            'A clip that stops short of the contig end would shorten it.'
        )
        min_clip = 1

    if exclude_outliers:
        logging.warning(
            '--exclude-outliers is deprecated and no longer excludes anything. '
            'Contigs with anomalous overhang coverage are reported in the '
            'stats report; exclude them with --exclude-contigs if you agree '
            'with the assessment.'
        )

    # Bound before the try so the finally clause can always close it.
    overhang_log_handle = None

    try:
        # Validate indexed files
        logging.info('Validating indexed input files...')
        is_valid, error_msg = validate_indexed_files(reference_fasta, bam_file)
        if not is_valid:
            raise click.ClickException(error_msg)

        # Set ref_idx as the reference FASTA with additional suffix .fai
        ref_idx = reference_fasta.parent / (reference_fasta.name + '.fai')

        # Handle default output filenames if not specified
        if not output_fasta and not dry_run:
            # Default to stdout (will be handled in writer logic)
            output_fasta = None

        if not stats_report:
            # Default stats to stderr (already handled in existing logic)
            stats_report = None

        # Validate output directories
        if output_fasta or stats_report:
            validate_output_directories(output_fasta, stats_report)

        logging.info('Reading reference genome index...')
        contig_dict = read_fai(ref_idx)
        logging.info(f'Loaded {len(contig_dict)} contigs from reference')

        # Parse excluded contigs from both string and file sources
        excluded_contig_set = combine_excluded_contigs(
            exclude_contigs, exclude_contigs_file, contig_dict
        )

        # Prepare motif patterns if specified
        motif_patterns = {}
        if count_motifs:
            motif_patterns = get_motif_regex(count_motifs, fuzzy_count)

        # Perform terminal motif screening if requested
        terminal_motif_counts = {}
        if screen_terminal_bases > 0 and motif_patterns:
            logging.info(
                f'Screening {screen_terminal_bases} terminal bases for motifs...'
            )
            terminal_motif_counts = count_terminal_motifs(
                reference_fasta, contig_dict, motif_patterns, screen_terminal_bases
            )
            total_screened = len(terminal_motif_counts)
            logging.info(
                f'Completed terminal motif screening for {total_screened} contigs'
            )

        # Open the per-overhang log, if requested, before streaming begins.
        if overhang_log:
            overhang_log.parent.mkdir(parents=True, exist_ok=True)
            overhang_log_handle = overhang_log.open('w', encoding='utf-8')
            overhang_log_handle.write('\t'.join(OVERHANG_LOG_HEADER) + '\n')
            logging.info(f'Writing per-overhang log to {overhang_log}')

        # Open indexed files
        logging.info('Opening indexed BAM and FASTA files...')
        with StreamingGenomeProcessor(reference_fasta, bam_file) as processor:
            bam_file_handle = pysam.AlignmentFile(str(bam_file), 'rb')

            # Stream contigs for extension analysis
            logging.info('Streaming contigs for extension analysis...')
            extensions_applied = {}
            excluded_contigs = []
            warnings = []
            motif_stats = {}
            post_motif_counts = {}
            all_stats = {}
            # Only populated when --html-report is requested.
            terminal_sequences = {}
            selected_reads = {}
            html_panels = []

            # Collect statistics for contigs that meet extension criteria
            extension_results = {}  # Store ExtensionResult objects for writing phase
            for contig_name, contig_stats in stream_contigs_for_extension(
                bam_file_handle,
                contig_dict,
                min_overhangs=min_overhangs,
                max_break=max_break,
                min_clip=min_clip,
                min_anchor=min_anchor,
                anchor_context=HTML_READ_CONTEXT if html_report else 0,
            ):
                all_stats[contig_name] = contig_stats
                logging.info(
                    f'{contig_name}: {contig_stats.left_count} left, '
                    f'{contig_stats.right_count} right overhang reads'
                )

                if overhang_log_handle is not None:
                    _write_overhang_log_rows(overhang_log_handle, contig_stats)

                # Check if this contig is explicitly excluded
                if contig_name in excluded_contig_set:
                    total_overhangs = contig_stats.left_count + contig_stats.right_count
                    logging.info(
                        f'Excluding contig "{contig_name}" (found {total_overhangs} overhangs: '
                        f'{contig_stats.left_count} left, {contig_stats.right_count} right)'
                    )
                    excluded_contigs.append(contig_name)
                    continue

                # Get the original sequence for this contig
                try:
                    original_sequence = processor.get_contig_sequence(contig_name)
                except KeyError:
                    logging.warning(
                        f'Contig {contig_name} not found in FASTA file, skipping'
                    )
                    continue

                if html_report:
                    # Keep only the terminal windows; the whole sequence is far
                    # too large to hold for every contig in an assembly.
                    window = min(HTML_CONTIG_CONTEXT, len(original_sequence))
                    terminal_sequences[contig_name] = (
                        original_sequence[:window].upper(),
                        original_sequence[-window:].upper(),
                    )

                # Process extension for this contig
                extension_result = process_single_contig_extension(
                    contig_name=contig_name,
                    contig_stats=contig_stats,
                    original_sequence=original_sequence,
                    min_extension=min_extension,
                    max_homopolymer=max_homopolymer,
                    motif_patterns=motif_patterns,
                    dry_run=dry_run,
                    terminal_length=screen_terminal_bases,
                )

                if html_report:
                    # Only the reads that will actually be rendered need their
                    # sequence kept. Pruning here bounds peak memory at roughly
                    # max_reads per end rather than every overhang in the
                    # assembly.
                    _prune_read_slices(contig_stats, max(1, html_max_reads))

                if extension_result:
                    # Store the complete ExtensionResult for later use
                    extension_results[contig_name] = extension_result
                    extensions_applied[contig_name] = extension_result.extension_info

                    if html_report:
                        # Record which read won at each end, so the alignment
                        # view can mark it among the candidates it beat.
                        info = extension_result.extension_info
                        selected_reads[contig_name] = {
                            'left': info.get('left_read_name')
                            if info.get('has_left_extension')
                            else None,
                            'right': info.get('right_read_name')
                            if info.get('has_right_extension')
                            else None,
                        }
                    warnings.extend(extension_result.warnings)
                    if extension_result.motif_counts:
                        motif_stats[contig_name] = extension_result.motif_counts
                    if extension_result.end_motif_counts:
                        post_motif_counts[contig_name] = (
                            extension_result.end_motif_counts
                        )

                    # Log successful extension(s)
                    ext_info = extension_result.extension_info

                    # Log both extensions if present
                    if ext_info.get('has_left_extension', False):
                        left_length = ext_info.get(
                            'left_overhang_length', ext_info.get('overhang_length', 0)
                        )
                        left_read = ext_info.get(
                            'left_read_name', ext_info.get('read_name', 'unknown')
                        )
                        if dry_run:
                            logging.info(
                                f'[DRY RUN] Would extend {contig_name} left end: '
                                f'+{left_length}bp from read {left_read}'
                            )
                        else:
                            logging.info(
                                f'Extended {contig_name} left end: '
                                f'+{left_length}bp from read {left_read}'
                            )

                    if ext_info.get('has_right_extension', False):
                        right_length = ext_info.get('right_overhang_length', 0)
                        right_read = ext_info.get('right_read_name', 'unknown')
                        if dry_run:
                            logging.info(
                                f'[DRY RUN] Would extend {contig_name} right end: '
                                f'+{right_length}bp from read {right_read}'
                            )
                        else:
                            logging.info(
                                f'Extended {contig_name} right end: '
                                f'+{right_length}bp from read {right_read}'
                            )

            if overhang_log_handle is not None:
                overhang_log_handle.close()
                overhang_log_handle = None

            # Calculate overall statistics if we have data
            if all_stats:
                overall_stats = calculate_overhang_statistics(all_stats)
                logging.info(f'Processed {len(all_stats)} contigs with overhangs')

                # Show the shape of the per-contig-end distribution. A long
                # right tail is the signature of a collapsed repeat or an
                # organellar contig attracting reads from across the assembly.
                end_counts = [s.left_count for s in all_stats.values()]
                end_counts += [s.right_count for s in all_stats.values()]
                chart = histogram(end_counts, label='overhang reads')
                if chart:
                    logging.info(
                        'Overhang reads per contig end '
                        f'({len(end_counts)} ends):\n{chart}'
                    )
            else:
                overall_stats = {'left': {}, 'right': {}}

            # Flag contigs whose overhang coverage stands out from their peers.
            # These are reported for review, never excluded automatically.
            anomalous = flag_anomalous_overhang_coverage(all_stats, outlier_threshold)
            flagged = sorted(
                set(anomalous['left_outliers']) | set(anomalous['right_outliers'])
            )
            if flagged:
                logging.warning(
                    f'{len(flagged)} contig end(s) have anomalous overhang '
                    'coverage, which can indicate a collapsed repeat, an rDNA '
                    'array or an organellar contig attracting reads from '
                    'elsewhere. Review the stats report and re-run with '
                    '--exclude-contigs if extension is not appropriate: '
                    + ', '.join(flagged)
                )

            if html_report:
                from ..report.html import render_contig_panels

                left_flags = set(anomalous['left_outliers'])
                right_flags = set(anomalous['right_outliers'])
                motif_list = sorted(motif_patterns) if motif_patterns else ()

                for contig_name, contig_stats in all_stats.items():
                    html_panels.extend(
                        render_contig_panels(
                            contig=contig_name,
                            stats=contig_stats,
                            terminal_sequences=terminal_sequences.get(
                                contig_name, ('', '')
                            ),
                            selected_reads=selected_reads.get(contig_name, {}),
                            flagged_left=contig_name in left_flags,
                            flagged_right=contig_name in right_flags,
                            motifs=motif_list,
                            max_reads=max(1, html_max_reads),
                            max_overhang=HTML_MAX_OVERHANG,
                            min_window=HTML_MIN_WINDOW,
                            window_cap=HTML_WINDOW_CAP,
                        )
                    )
                    # The rendered blocks are all that is needed from here on.
                    for oh in (
                        contig_stats.left_overhangs + contig_stats.right_overhangs
                    ):
                        oh.read_seq = ''

            # Write extended sequences
            if not dry_run:
                if output_fasta:
                    logging.info(f'Writing extended sequences to {output_fasta}...')
                else:
                    logging.info('Writing extended sequences to stdout...')
                # Write in a single pass over the reference index, so output
                # contig order matches input order. Writing all extended
                # contigs first and appending the rest afterwards would shunt
                # every excluded, unsupported or failed contig to the end of
                # the file, breaking any downstream tool that assumes
                # input/output correspondence.
                with BufferedContigWriter(output_fasta) as writer:
                    missing = 0
                    for contig_name in contig_dict:
                        extension_result = extension_results.get(contig_name)
                        if extension_result is not None:
                            writer.write_contig(
                                contig_name, extension_result.extended_sequence
                            )
                            logging.debug(f'Wrote extended sequence for {contig_name}')
                            continue

                        try:
                            writer.write_contig(
                                contig_name,
                                processor.get_contig_sequence(contig_name),
                            )
                        except KeyError:
                            # Present in the .fai but not the FASTA itself.
                            missing += 1
                            logging.warning(
                                f'Contig {contig_name} is in the index but not the '
                                'FASTA; omitted from output'
                            )

                    if missing:
                        logging.warning(
                            f'{missing} indexed contig(s) were missing from the '
                            'FASTA and are absent from the output'
                        )

            bam_file_handle.close()

        # Generate report
        report_content = generate_extension_report(
            all_stats,
            extensions_applied,
            anomalous,
            overall_stats,
            excluded_contigs,
            warnings,
            motif_stats,
            terminal_motif_counts,
            dry_run,
            post_motif_counts=post_motif_counts,
            screen_terminal_bases=screen_terminal_bases,
            total_contigs=len(contig_dict),
        )

        # Write outputs
        if stats_report:
            if str(stats_report) == '-':
                # Write to stdout when explicitly requested with '-'
                logging.info('Writing statistics report to stdout')
                print(report_content)
            else:
                logging.info(f'Writing statistics report to {stats_report}')
                with open(stats_report, 'w') as f:
                    f.write(report_content)
        else:
            # Default: write report to stderr if no file specified
            logging.info('Writing statistics report to stderr')
            print(report_content, file=sys.stderr)

        if html_report:
            from .._version import __version__
            from ..report.html import render_html_report

            logging.info(f'Writing HTML report to {html_report}...')
            html_report.parent.mkdir(parents=True, exist_ok=True)
            html_report.write_text(
                render_html_report(
                    stats_dict=all_stats,
                    extensions_applied=extensions_applied,
                    anomalous=anomalous,
                    excluded_contigs=excluded_contigs,
                    warnings=warnings,
                    panels=html_panels,
                    total_contigs=len(contig_dict),
                    dry_run=dry_run,
                    motifs=sorted(motif_patterns) if motif_patterns else (),
                    max_reads=max(1, html_max_reads),
                    max_overhang=HTML_MAX_OVERHANG,
                    version=__version__,
                    # Verbatim, so the report can be reproduced from itself.
                    command=shlex.join(sys.argv),
                    generated=datetime.now()
                    .astimezone()
                    .strftime('%Y-%m-%d %H:%M:%S %Z'),
                ),
                encoding='utf-8',
            )
            logging.info(f'HTML report written to {html_report}')

        # Summary
        if dry_run:
            logging.info(
                f'Dry-run analysis complete: {len(extensions_applied)} contigs would be extended'
            )
        else:
            logging.info(
                f'Extension complete: {len(extensions_applied)} contigs extended'
            )

            # Add polishing reminder for actual extensions
            if extensions_applied and not dry_run:
                logging.info(
                    'IMPORTANT: Extended contigs should be polished with appropriate tools '
                    '(e.g., Medaka for ONT data, and Pypolca for Illumina data) to improve accuracy before downstream analysis.'
                )

        # Motif analysis summary
        if motif_patterns and motif_stats:
            total_motifs_found = sum(
                sum(counts.values()) for counts in motif_stats.values()
            )
            contigs_with_motifs = len(
                [
                    contig
                    for contig, counts in motif_stats.items()
                    if any(count > 0 for count in counts.values())
                ]
            )
            if total_motifs_found > 0:
                logging.info(
                    f'Motif analysis: found {total_motifs_found} motif matches '
                    f'in {contigs_with_motifs} extended contigs'
                )

        if excluded_contigs:
            logging.info(f'Excluded {len(excluded_contigs)} outlier contigs')
        if warnings:
            logging.info(f'Generated {len(warnings)} warnings')

    except click.ClickException:
        # Already carries a user-facing message; do not wrap it again.
        raise
    except Exception as e:
        logging.error(f'Error during extend operation: {e}')
        raise click.ClickException(str(e)) from e
    finally:
        # Close the overhang log even if the run failed partway through, so
        # whatever was written is readable.
        if overhang_log_handle is not None:
            overhang_log_handle.close()
