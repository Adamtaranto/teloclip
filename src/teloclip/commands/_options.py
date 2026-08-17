"""
Argument parsing and validation helpers shared by the ``extend`` sub-command.

These turn what the user typed into what the pipeline needs: output paths
checked before any work begins, exclusion lists gathered from repeated options
and files, and motif strings compiled to regexes. They are kept out of the
command module so that the command reads as a sequence of steps rather than a
mixture of parsing and orchestration.

Failures here raise ``click`` exceptions: a bad path or an unreadable exclusion
file is a usage error, and the user should see a message rather than a
traceback.
"""

import logging
from pathlib import Path
import re
from typing import Dict, List

import click

from ..core.motifs import make_fuzzy_motif_regex, make_motif_regex
from ..core.seqops import revComp


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
