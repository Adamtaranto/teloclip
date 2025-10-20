"""
Extend sub-command implementation.

This module implements the 'teloclip extend' command for automatically extending
draft contigs using overhang analysis from soft-clipped alignments.

This version is optimized for large genomes using streaming I/O and indexed access
to avoid loading entire genomes into memory.
"""

import logging
from pathlib import Path
import re
import sys
from typing import Dict, List

import click
import pysam

from ..analysis import ContigStats, calculate_overhang_statistics
from ..motifs import make_fuzzy_motif_regex, make_motif_regex
from ..seqops import read_fai, revComp
from ..streaming_analysis import (
    process_single_contig_extension,
    stream_contigs_for_extension,
)
from ..streaming_io import (
    BufferedContigWriter,
    StreamingGenomeProcessor,
    copy_unmodified_contigs,
    validate_indexed_files,
)


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


def generate_extension_report(
    stats_dict: Dict[str, ContigStats],
    extensions_applied: Dict[str, dict],
    outliers: Dict[str, List[str]],
    overall_stats: Dict[str, Dict[str, float]],
    excluded_contigs: List[str],
    warnings: List[str],
    motif_stats: Dict[str, Dict[str, int]] = None,
    dry_run: bool = False,
) -> str:
    """
    Generate a comprehensive statistics report for contig extension analysis.

    Parameters
    ----------
    stats_dict : Dict[str, ContigStats]
        Dictionary mapping contig names to their statistics.
    extensions_applied : Dict[str, dict]
        Dictionary of extensions that were applied to contigs.
    outliers : Dict[str, List[str]]
        Dictionary of outlier contigs by category.
    overall_stats : Dict[str, Dict[str, float]]
        Overall statistics across all contigs.
    excluded_contigs : List[str]
        List of contig names that were excluded from analysis.
    warnings : List[str]
        List of warning messages generated during analysis.
    motif_stats : Dict[str, Dict[str, int]], optional
        Statistics about motif occurrences. Default is None.
    dry_run : bool, optional
        Whether this is a dry run (no actual extensions applied). Default is False.

    Returns
    -------
    str
        Formatted statistics report as a multi-line string.
    """
    report_lines = []

    if dry_run:
        report_lines.append('# Teloclip Extend Statistics Report (DRY RUN)')
    else:
        report_lines.append('# Teloclip Extend Statistics Report')
    report_lines.append('=' * 50)
    report_lines.append('')

    # Overall statistics
    report_lines.append('## Overall Overhang Statistics')
    for category, stats in overall_stats.items():
        report_lines.append(f'\n### {category.title()} Overhangs')
        for metric, value in stats.items():
            report_lines.append(f'  {metric}: {value:.2f}')

    report_lines.append('')

    # Extensions applied
    if dry_run:
        report_lines.append('## Extensions That Would Be Applied')
        report_lines.append(
            f'Total contigs that would be extended: {len(extensions_applied)}'
        )
    else:
        report_lines.append('## Extensions Applied')
        report_lines.append(f'Total contigs extended: {len(extensions_applied)}')
    report_lines.append('')

    for contig_name, ext_info in extensions_applied.items():
        report_lines.append(f'### {contig_name}')
        report_lines.append(
            f'  Direction: {"Left" if ext_info["is_left"] else "Right"}'
        )
        report_lines.append(f'  Original length: {ext_info["original_length"]:,}')
        report_lines.append(f'  Extension length: {ext_info["overhang_length"]}')
        report_lines.append(f'  Final length: {ext_info["final_length"]:,}')
        report_lines.append(f'  Source read: {ext_info["read_name"]}')
        if ext_info['trim_length'] > 0:
            report_lines.append(f'  Bases trimmed: {ext_info["trim_length"]}')
        report_lines.append('')

    # Motif analysis results
    if motif_stats:
        report_lines.append('## Motif Analysis Results')
        report_lines.append('')
        for contig_name, motif_counts in motif_stats.items():
            if any(count > 0 for count in motif_counts.values()):
                report_lines.append(f'### {contig_name}')
                for motif_name, count in motif_counts.items():
                    report_lines.append(f'  {motif_name}: {count} matches')
                report_lines.append('')

    # Outliers detected
    if any(outliers.values()):
        report_lines.append('## Outlier Contigs Detected')
        if outliers['left_outliers']:
            report_lines.append(
                f'Left outliers: {", ".join(outliers["left_outliers"])}'
            )
        if outliers['right_outliers']:
            report_lines.append(
                f'Right outliers: {", ".join(outliers["right_outliers"])}'
            )
        report_lines.append('')

    # Excluded contigs
    if excluded_contigs:
        report_lines.append('## Excluded Contigs')
        for contig in excluded_contigs:
            report_lines.append(f'  {contig}')
        report_lines.append('')

    # Warnings
    if warnings:
        report_lines.append('## Warnings')
        for warning in warnings:
            report_lines.append(f'  - {warning}')
        report_lines.append('')

    # Per-contig summary
    report_lines.append('## Per-Contig Overhang Summary')
    report_lines.append(
        'Contig\tLength\tLeft_Count\tRight_Count\tLeft_Total\tRight_Total'
    )

    for contig_name, contig_stats in stats_dict.items():
        report_lines.append(
            f'{contig_name}\t{contig_stats.contig_length}\t'
            f'{contig_stats.left_count}\t{contig_stats.right_count}\t'
            f'{contig_stats.left_total_length}\t{contig_stats.right_total_length}'
        )

    return '\n'.join(report_lines)


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
                pattern = make_fuzzy_motif_regex(motif)
                pattern_name = f'{motif}_fuzzy'
                logging.debug(f'Created fuzzy pattern for {motif}: {pattern}')
            else:
                pattern = make_motif_regex(motif)
                pattern_name = motif
                logging.debug(f'Created exact pattern for {motif}: {pattern}')

            motif_patterns[pattern_name] = pattern

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
    '--exclude-outliers', is_flag=True, help='Exclude outlier contigs from extension'
)
@click.option(
    '--outlier-threshold',
    type=float,
    default=2.0,
    help='Z-score threshold for outlier detection (default: 2.0)',
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
    default=50,
    help='Maximum homopolymer run length allowed (default: 50)',
)
@click.option(
    '--min-extension',
    type=int,
    default=1,
    help='Minimum overhang length for extension (default: 1)',
)
@click.option(
    '--max-break',
    type=int,
    default=10,
    help='Maximum gap allowed between alignment and contig end (default: 10)',
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
    '--log-level',
    default='INFO',
    type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR'], case_sensitive=False),
    help='Logging level (default: INFO).',
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
    max_break,
    min_anchor,
    dry_run,
    count_motifs,
    fuzzy_count,
    log_level,
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
        Exclude outlier overhangs from analysis.
    outlier_threshold : float
        Z-score threshold for outlier detection.
    min_overhangs : int
        Minimum number of overhangs required for extension.
    max_homopolymer : int
        Maximum homopolymer length allowed in extensions.
    min_extension : int
        Minimum extension length to report.
    max_break : int
        Maximum distance from contig end to search for overhangs.
    min_anchor : int
        Minimum anchor length in aligned portion.
    dry_run : bool
        Perform analysis without writing output files.
    count_motifs : bool
        Count telomeric motifs in extensions.
    fuzzy_count : int
        Number of mismatches allowed in motif counting.
    log_level : str
        Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
    """
    from ..logs import init_logging

    # Initialize logging for this command
    init_logging(log_level)

    ctx.ensure_object(dict)

    try:
        # Validate indexed files
        logging.info('Validating indexed input files...')
        is_valid, error_msg = validate_indexed_files(reference_fasta, bam_file)
        if not is_valid:
            raise click.ClickException(error_msg)

        # Set ref_idx as the reference FASTA with additional suffix .fai
        ref_idx = reference_fasta.parent / (reference_fasta.name + '.fai')

        # Validate output directories
        if output_fasta or stats_report:
            validate_output_directories(output_fasta, stats_report)

        logging.info('Reading reference genome index...')
        contig_dict = read_fai(ref_idx)
        logging.info(f'Loaded {len(contig_dict)} contigs from reference')

        # Prepare motif patterns if specified
        motif_patterns = {}
        if count_motifs:
            motif_patterns = get_motif_regex(count_motifs, fuzzy_count)

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
            all_stats = {}

            # Collect statistics for contigs that meet extension criteria
            for contig_name, contig_stats in stream_contigs_for_extension(
                bam_file_handle,
                contig_dict,
                min_overhangs=min_overhangs,
                max_break=max_break,
                min_clip=1,  # Use default minimum clip
                min_anchor=min_anchor,
                exclude_outliers=exclude_outliers,
                outlier_threshold=outlier_threshold,
            ):
                all_stats[contig_name] = contig_stats
                logging.debug(f'Processing contig {contig_name} for extension...')

                # Get the original sequence for this contig
                try:
                    original_sequence = processor.get_contig_sequence(contig_name)
                except KeyError:
                    logging.warning(
                        f'Contig {contig_name} not found in FASTA file, skipping'
                    )
                    continue

                # Process extension for this contig
                extension_result = process_single_contig_extension(
                    contig_name=contig_name,
                    contig_stats=contig_stats,
                    original_sequence=original_sequence,
                    min_extension=min_extension,
                    max_homopolymer=max_homopolymer,
                    motif_patterns=motif_patterns,
                    dry_run=dry_run,
                )

                if extension_result:
                    extensions_applied[contig_name] = extension_result.extension_info
                    warnings.extend(extension_result.warnings)
                    if extension_result.motif_counts:
                        motif_stats[contig_name] = extension_result.motif_counts

                    # Log successful extension
                    ext_info = extension_result.extension_info
                    end_name = 'left' if ext_info['is_left'] else 'right'
                    if dry_run:
                        logging.info(
                            f'[DRY RUN] Would extend {contig_name} {end_name} end: '
                            f'+{ext_info["overhang_length"]}bp from read {ext_info["read_name"]}'
                        )
                    else:
                        logging.info(
                            f'Extended {contig_name} {end_name} end: '
                            f'+{ext_info["overhang_length"]}bp from read {ext_info["read_name"]}'
                        )

            # Calculate overall statistics if we have data
            if all_stats:
                overall_stats = calculate_overhang_statistics(all_stats)
                logging.info(f'Processed {len(all_stats)} contigs with overhangs')
            else:
                overall_stats = {'left': {}, 'right': {}}

            # Write extended sequences
            if not dry_run:
                logging.info('Writing extended sequences...')
                with BufferedContigWriter(output_fasta) as writer:
                    # Write extended contigs
                    extended_contig_names = set()
                    for contig_name, extension_result in [
                        (name, er)
                        for name, er in [
                            (
                                n,
                                process_single_contig_extension(
                                    n,
                                    all_stats[n],
                                    processor.get_contig_sequence(n),
                                    min_extension,
                                    max_homopolymer,
                                    motif_patterns,
                                    dry_run,
                                ),
                            )
                            for n in extensions_applied.keys()
                        ]
                        if er is not None
                    ]:
                        writer.write_contig(
                            contig_name, extension_result.extended_sequence
                        )
                        extended_contig_names.add(contig_name)

                    # Copy unmodified contigs
                    copy_unmodified_contigs(
                        processor, writer, extended_contig_names, contig_dict
                    )

            bam_file_handle.close()

        # Generate report
        report_content = generate_extension_report(
            all_stats,
            extensions_applied,
            {
                'left_outliers': [],
                'right_outliers': [],
            },  # Outlier detection handled in streaming
            overall_stats,
            excluded_contigs,
            warnings,
            motif_stats,
            dry_run,
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

        # Summary
        if dry_run:
            logging.info(
                f'Dry-run analysis complete: {len(extensions_applied)} contigs would be extended'
            )
        else:
            logging.info(
                f'Extension complete: {len(extensions_applied)} contigs extended'
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

    except Exception as e:
        logging.error(f'Error during extend operation: {e}')
        raise click.ClickException(str(e)) from e
