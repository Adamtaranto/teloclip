"""
Extend sub-command implementation.

This module implements the 'teloclip extend' command for automatically extending
draft contigs using overhang analysis from soft-clipped alignments.
"""

# TODO: Use biopython seqIO for fasta reading/writing so that we can handle large genomes without loading everything into memory
## This will require changing or removing fasta2dict and writefasta functions to use generators. Check for other code that assumes fasta2dict returns full dict in memory.
# TODO: Skip contigs not present in reference fasta when reading fai index and warn user
# TODO: Write output fasta in streaming fashion to stdout if no output file specified
# TODO: Option to read sam from stdin if sam_file is '-'
# TODO: Do not report extensions for contigs not present in reference fasta
# TODO: Do not report extensions complete in dry-run mode
# TODO: Fix error error on extend.py:397 Error during extend operation: can only concatenate tuple (not "str") to tuple. Error: can only concatenate tuple (not "str") to tuple
# TODO: Default stats-report to stderr if not specified
# TODO: Report count of motif matches in each extended region using re.findall

import logging
from pathlib import Path
import sys
from typing import Dict, Iterator, List

import click

from ..analysis import (
    ContigStats,
    calculate_overhang_statistics,
    collect_overhang_stats,
    detect_homopolymer_runs,
    identify_outlier_contigs,
    select_best_overhang,
)
from ..extension import apply_contig_extension
from ..logs import init_logging
from ..bio_io import load_fasta_sequences, write_fasta_sequences, validate_fasta_against_fai
from ..seqops import read_fai


def setup_logger(level):
    """Setup logger with the specified level."""
    init_logging()
    logger = logging.getLogger()
    logger.setLevel(level)
    return logger


def validate_input_files(sam_file: Path, reference_fasta: Path, ref_idx: Path) -> None:
    """Validate that required input files exist and are readable."""
    if not sam_file.exists():
        raise click.ClickException(f'SAM file not found: {sam_file}')

    if not reference_fasta.exists():
        raise click.ClickException(f'Reference FASTA not found: {reference_fasta}')

    if not ref_idx.exists():
        raise click.ClickException(f'Reference index not found: {ref_idx}')


def validate_output_directories(output_fasta: Path, stats_report: Path) -> None:
    """Validate that output directories are writable."""
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


def read_sam_lines(sam_file: Path) -> Iterator[str]:
    """Read SAM file lines, handling both files and stdin."""
    if str(sam_file) == '-':
        for line in sys.stdin:
            yield line
    else:
        with open(sam_file, 'r') as f:
            for line in f:
                yield line


def generate_extension_report(
    stats_dict: Dict[str, ContigStats],
    extensions_applied: Dict[str, dict],
    outliers: Dict[str, List[str]],
    overall_stats: Dict[str, Dict[str, float]],
    excluded_contigs: List[str],
    warnings: List[str],
) -> str:
    """Generate a comprehensive statistics report."""
    report_lines = []

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


@click.command()
@click.argument('sam_file', type=click.Path(exists=False, path_type=Path))
@click.argument('reference_fasta', type=click.Path(exists=True, path_type=Path))
@click.option(
    '--ref-idx',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Path to fai index for reference fasta',
)
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
    '--dry-run', is_flag=True, help='Report extensions without modifying sequences'
)
@click.pass_context
def extend(
    ctx,
    sam_file,
    reference_fasta,
    ref_idx,
    output_fasta,
    stats_report,
    exclude_outliers,
    outlier_threshold,
    min_overhangs,
    max_homopolymer,
    min_extension,
    dry_run,
):
    """
    Extend contigs using overhang analysis from soft-clipped alignments.

    This command analyzes soft-clipped alignments to identify overhanging sequences
    that extend beyond contig ends, then automatically extends contigs using the
    longest suitable overhangs.

    SAM_FILE can be a file path or '-' to read from stdin.
    REFERENCE_FASTA should be the original reference used for alignment.
    """
    logger = setup_logger(ctx.obj.get('log_level', 'INFO'))

    try:
        # Validate input files
        validate_input_files(sam_file, reference_fasta, ref_idx)

        # Validate output directories
        if output_fasta or stats_report:
            validate_output_directories(output_fasta, stats_report)

        logger.info('Reading reference genome index...')
        contig_dict = read_fai(ref_idx)
        logger.info(f'Loaded {len(contig_dict)} contigs from reference')

        logger.info('Reading reference sequences...')
        reference_seqs = load_fasta_sequences(reference_fasta)

        logger.info('Collecting overhang statistics...')
        sam_lines = read_sam_lines(sam_file)
        stats_dict = collect_overhang_stats(sam_lines, contig_dict)

        # Calculate overall statistics
        overall_stats = calculate_overhang_statistics(stats_dict)
        logger.info(
            f'Collected overhangs from {sum(s.left_count + s.right_count for s in stats_dict.values())} alignments'
        )

        # Identify outliers
        outliers = identify_outlier_contigs(stats_dict, outlier_threshold)
        logger.info(
            f'Identified {len(outliers["left_outliers"])} left outliers and {len(outliers["right_outliers"])} right outliers'
        )

        # Track extensions and warnings
        extensions_applied = {}
        excluded_contigs = []
        warnings = []

        # Process each contig for potential extension
        for contig_name, contig_stats in stats_dict.items():
            # Skip if excluding outliers
            if exclude_outliers and (
                contig_name in outliers['left_outliers']
                or contig_name in outliers['right_outliers']
            ):
                excluded_contigs.append(contig_name)
                logger.debug(f'Skipping outlier contig: {contig_name}')
                continue

            # Check each end for extension opportunities
            for is_left in [True, False]:
                overhangs = (
                    contig_stats.left_overhangs
                    if is_left
                    else contig_stats.right_overhangs
                )
                end_name = 'left' if is_left else 'right'

                # Skip if insufficient overhangs
                if len(overhangs) < min_overhangs:
                    logger.debug(
                        f'Insufficient {end_name} overhangs for {contig_name}: {len(overhangs)} < {min_overhangs}'
                    )
                    continue

                # Select best overhang
                best_overhang = select_best_overhang(
                    overhangs, min_extension, max_homopolymer
                )

                if best_overhang is None:
                    logger.debug(
                        f'No suitable {end_name} overhang found for {contig_name}'
                    )
                    continue

                # Check for homopolymer runs
                homo_runs = detect_homopolymer_runs(
                    best_overhang.sequence, max_homopolymer
                )
                if homo_runs:
                    warning_msg = (
                        f'Homopolymer run detected in {contig_name} {end_name} extension: '
                        f'{homo_runs[0][0]}x{homo_runs[0][3]} at position {homo_runs[0][1]}'
                    )
                    warnings.append(warning_msg)
                    logger.warning(warning_msg)

                # Apply extension (only one per contig for now)
                if contig_name not in extensions_applied:
                    if not dry_run:
                        try:
                            # Extract sequence from tuple (header, seq)
                            original_seq = reference_seqs[contig_name][1]
                            extended_seq, ext_info = apply_contig_extension(
                                original_seq, best_overhang, contig_stats.contig_length
                            )
                            # Update sequence in tuple (header, seq)
                            header = reference_seqs[contig_name][0]
                            reference_seqs[contig_name] = (header, extended_seq)
                            extensions_applied[contig_name] = ext_info

                            logger.info(
                                f'Extended {contig_name} {end_name} end: '
                                f'+{best_overhang.length}bp from read {best_overhang.read_name}'
                            )
                        except ValueError as e:
                            error_msg = f'Extension failed for {contig_name}: {e}'
                            warnings.append(error_msg)
                            logger.error(error_msg)
                    else:
                        # Dry run mode
                        ext_info = {
                            'overhang_length': best_overhang.length,
                            'read_name': best_overhang.read_name,
                            'is_left': is_left,
                            'original_length': contig_stats.contig_length,
                            'final_length': contig_stats.contig_length
                            + best_overhang.length,
                            'trim_length': 0,  # Simplified for dry run
                        }
                        extensions_applied[contig_name] = ext_info

                        logger.info(
                            f'[DRY RUN] Would extend {contig_name} {end_name} end: '
                            f'+{best_overhang.length}bp from read {best_overhang.read_name}'
                        )

        # Generate report
        report_content = generate_extension_report(
            stats_dict,
            extensions_applied,
            outliers,
            overall_stats,
            excluded_contigs,
            warnings,
        )

        # Write outputs
        if stats_report:
            logger.info(f'Writing statistics report to {stats_report}')
            with open(stats_report, 'w') as f:
                f.write(report_content)
        else:
            # Print report to stdout if no file specified
            print(report_content)

        if output_fasta and not dry_run:
            logger.info(f'Writing extended sequences to {output_fasta}')
            write_fasta_sequences(reference_seqs, output_fasta)

        # Summary
        logger.info(f'Extension complete: {len(extensions_applied)} contigs extended')
        if excluded_contigs:
            logger.info(f'Excluded {len(excluded_contigs)} outlier contigs')
        if warnings:
            logger.info(f'Generated {len(warnings)} warnings')

    except Exception as e:
        logger.error(f'Error during extend operation: {e}')
        raise click.ClickException(str(e)) from e
