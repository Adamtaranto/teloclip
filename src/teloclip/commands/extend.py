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
from typing import Dict

import click
import pyfaidx
import pysam

from ..core.analysis import (
    calculate_overhang_statistics,
    flag_anomalous_overhang_coverage,
)
from ..core.seqops import read_fai
from ..core.streaming_analysis import (
    process_single_contig_extension,
    stream_contigs_for_extension,
)
from ..io.streaming import (
    BufferedContigWriter,
    StreamingGenomeProcessor,
    validate_indexed_files,
)
from ..report.extend_report import (
    HTML_CONTIG_CONTEXT,
    HTML_MAX_OVERHANG,
    HTML_MIN_WINDOW,
    HTML_READ_CONTEXT,
    HTML_WINDOW_CAP,
    OVERHANG_LOG_HEADER,
    generate_extension_report,
    prune_read_slices,
    write_overhang_log_rows,
)
from ..report.text import histogram
from ._options import (
    combine_excluded_contigs,
    get_motif_regex,
    read_excluded_contigs_file,
    validate_output_directories,
)

# Re-exported rather than merely used. The argument-parsing helpers and the
# report builder were defined in this module before it was split, and tests
# and external callers still reach for them at teloclip.commands.extend.
__all__ = [
    'combine_excluded_contigs',
    'count_terminal_motifs',
    'extend',
    'generate_extension_report',
    'get_motif_regex',
    'read_excluded_contigs_file',
    'validate_output_directories',
]


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
                    write_overhang_log_rows(overhang_log_handle, contig_stats)

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
                    prune_read_slices(contig_stats, max(1, html_max_reads))

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

            # Length is scored separately from depth. Reported at info rather
            # than warning level: unusually long overhangs on their own often
            # mean a genuine long telomere the assembly stopped short of, which
            # is the case extension exists to handle. It is the overlap with
            # the depth flags above that suggests a collapsed array.
            length_flagged = sorted(
                set(anomalous['left_length_outliers'])
                | set(anomalous['right_length_outliers'])
            )
            if length_flagged:
                also_deep = sorted(set(length_flagged) & set(flagged))
                logging.info(
                    f'{len(length_flagged)} contig end(s) have unusually long '
                    'overhangs relative to the rest of the assembly: '
                    + ', '.join(length_flagged)
                )
                if also_deep:
                    logging.warning(
                        'These contigs are anomalous in both overhang depth and '
                        'overhang length, the pattern typical of a collapsed '
                        'repeat or rDNA array at the contig end: '
                        + ', '.join(also_deep)
                    )

            if html_report:
                from ..report.panels import render_contig_panels

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
