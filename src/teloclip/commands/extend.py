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
import pyfaidx
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
    terminal_motif_counts: Dict[str, Dict[str, Dict[str, int]]] = None,
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
    terminal_motif_counts : Dict[str, Dict[str, Dict[str, int]]], optional
        Pre-extension terminal motif counts. Default is None.
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
        report_lines.append(f'  Original length: {ext_info["original_length"]:,}')
        report_lines.append(f'  Final length: {ext_info["final_length"]:,}')

        # Report left extension if present
        if ext_info.get('has_left_extension', False):
            left_length = ext_info.get('left_overhang_length', 0)
            left_read = ext_info.get('left_read_name', 'unknown')
            left_trim = ext_info.get('left_trim_length', 0)
            report_lines.append(
                f'  Left extension: +{left_length}bp from read {left_read}'
            )
            if left_trim > 0:
                report_lines.append(f'    Left bases trimmed: {left_trim}')

        # Report right extension if present
        if ext_info.get('has_right_extension', False):
            right_length = ext_info.get('right_overhang_length', 0)
            right_read = ext_info.get('right_read_name', 'unknown')
            right_trim = ext_info.get('right_trim_length', 0)
            report_lines.append(
                f'  Right extension: +{right_length}bp from read {right_read}'
            )
            if right_trim > 0:
                report_lines.append(f'    Right bases trimmed: {right_trim}')

        # Fallback for backward compatibility (single extension)
        if not ext_info.get('has_left_extension', False) and not ext_info.get(
            'has_right_extension', False
        ):
            direction = 'Left' if ext_info.get('is_left', False) else 'Right'
            report_lines.append(f'  Direction: {direction}')
            report_lines.append(
                f'  Extension length: {ext_info.get("overhang_length", 0)}'
            )
            report_lines.append(
                f'  Source read: {ext_info.get("read_name", "unknown")}'
            )
            if ext_info.get('trim_length', 0) > 0:
                report_lines.append(f'  Bases trimmed: {ext_info["trim_length"]}')

        report_lines.append('')

    # Terminal motif screening results
    if terminal_motif_counts:
        report_lines.append('## Terminal Region Motif Analysis (Pre-Extension)')
        report_lines.append('')
        for contig_name, terminal_counts in terminal_motif_counts.items():
            left_total = sum(terminal_counts['left'].values())
            right_total = sum(terminal_counts['right'].values())
            if left_total > 0 or right_total > 0:
                report_lines.append(f'### {contig_name}')
                report_lines.append(f'  Left terminal: {left_total} total motifs')
                for motif_name, count in terminal_counts['left'].items():
                    if count > 0:
                        report_lines.append(f'    {motif_name}: {count}')
                report_lines.append(f'  Right terminal: {right_total} total motifs')
                for motif_name, count in terminal_counts['right'].items():
                    if count > 0:
                        report_lines.append(f'    {motif_name}: {count}')
                report_lines.append('')

    # Extension motif analysis results
    if motif_stats:
        report_lines.append('## Extension Region Motif Analysis (Post-Extension)')
        report_lines.append('')
        for contig_name, motif_counts in motif_stats.items():
            if any(count > 0 for count in motif_counts.values()):
                report_lines.append(f'### {contig_name}')
                for motif_name, count in motif_counts.items():
                    if count > 0:
                        report_lines.append(f'  {motif_name}: {count} matches')

                        # Calculate gain if terminal counts available
                        if (
                            terminal_motif_counts
                            and contig_name in terminal_motif_counts
                        ):
                            # Determine which terminal was extended
                            extension_info = extensions_applied.get(contig_name, {})
                            if extension_info:
                                is_left = extension_info.get('is_left', False)
                                terminal_side = 'left' if is_left else 'right'
                                pre_count = terminal_motif_counts[contig_name][
                                    terminal_side
                                ].get(motif_name, 0)
                                gain = count - pre_count
                                if gain != 0:
                                    report_lines.append(
                                        f'    Gain from extension: {gain:+d}'
                                    )
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

                # Initialize counts for this contig
                terminal_counts[contig_name] = {
                    'left': {},
                    'right': {},
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


def parse_excluded_contigs(
    exclude_contigs_str: str, contig_dict: Dict[str, int]
) -> set:
    """
    Parse and validate excluded contig names.

    Parameters
    ----------
    exclude_contigs_str : str
        Comma-delimited string of contig names to exclude.
    contig_dict : Dict[str, int]
        Dictionary mapping contig names to their lengths.

    Returns
    -------
    set
        Set of valid contig names to exclude.
    """
    excluded_set = set()

    if not exclude_contigs_str:
        return excluded_set

    # Parse comma-delimited contig names
    raw_contigs = [
        contig.strip() for contig in exclude_contigs_str.split(',') if contig.strip()
    ]

    if not raw_contigs:
        return excluded_set

    logging.info(f'Processing exclusion list: {", ".join(raw_contigs)}')

    # Validate each contig name
    for contig_name in raw_contigs:
        if contig_name in contig_dict:
            excluded_set.add(contig_name)
            logging.info(f'Contig "{contig_name}" will be excluded from extension')
        else:
            logging.warning(
                f'Excluded contig "{contig_name}" not found in reference FASTA index'
            )

    if excluded_set:
        logging.info(f'Total contigs excluded: {len(excluded_set)}')
    else:
        logging.warning('No valid contigs found in exclusion list')

    return excluded_set


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
        Combined set of valid contig names to exclude.
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

    # Validate each contig name
    excluded_set = set()
    for contig_name in unique_names:
        if contig_name in contig_dict:
            excluded_set.add(contig_name)
            logging.info(f'Contig "{contig_name}" will be excluded from extension')
        else:
            logging.warning(
                f'Excluded contig "{contig_name}" not found in reference FASTA index'
            )

    if excluded_set:
        logging.info(f'Total unique contigs excluded: {len(excluded_set)}')
    else:
        logging.warning('No valid contigs found in exclusion sources')

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
    default=500,
    help='Maximum homopolymer run length allowed (default: 500)',
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
    prefix,
    screen_terminal_bases,
    exclude_contigs,
    exclude_contigs_file,
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
            extension_results = {}  # Store ExtensionResult objects for writing phase
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
                    # Store the complete ExtensionResult for later use
                    extension_results[contig_name] = extension_result
                    extensions_applied[contig_name] = extension_result.extension_info
                    warnings.extend(extension_result.warnings)
                    if extension_result.motif_counts:
                        motif_stats[contig_name] = extension_result.motif_counts

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

            # Calculate overall statistics if we have data
            if all_stats:
                overall_stats = calculate_overhang_statistics(all_stats)
                logging.info(f'Processed {len(all_stats)} contigs with overhangs')
            else:
                overall_stats = {'left': {}, 'right': {}}

            # Write extended sequences
            if not dry_run:
                if output_fasta:
                    logging.info(f'Writing extended sequences to {output_fasta}...')
                else:
                    logging.info('Writing extended sequences to stdout...')
                with BufferedContigWriter(output_fasta) as writer:
                    # Write extended contigs using stored results
                    extended_contig_names = set()
                    for contig_name, extension_result in extension_results.items():
                        writer.write_contig(
                            contig_name, extension_result.extended_sequence
                        )
                        extended_contig_names.add(contig_name)
                        logging.debug(f'Wrote extended sequence for {contig_name}')

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
            terminal_motif_counts,
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

    except Exception as e:
        logging.error(f'Error during extend operation: {e}')
        raise click.ClickException(str(e)) from e
