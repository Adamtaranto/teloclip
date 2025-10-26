"""
Extract sub-command for teloclip CLI.
"""

import logging
from pathlib import Path
import sys
from typing import Dict

import click

from ..extract_io import ExtractionStats
from ..logs import init_logging
from ..motifs import make_fuzzy_motif_regex, make_motif_regex
from ..samops import (
    EnhancedStreamingSamFilter,
    enhanced_streaming_split_by_contig,
)
from ..seqops import read_fai, revComp


@click.command(
    'extract',
    help='Extract overhanging reads for each end of each reference contig. Reads are always written to output files.',
)
@click.argument('samfile', type=click.File('r'), default=sys.stdin)
@click.option(
    '--ref-idx',
    required=True,
    type=click.Path(exists=True),
    help='Path to fai index for reference fasta. Index fasta using `samtools faidx FASTA`',
)
@click.option(
    '--prefix', type=str, help='Use this prefix for output files. Default: None.'
)
@click.option(
    '--extract-dir',
    type=click.Path(),
    help='Write extracted reads to this directory. Default: cwd.',
)
@click.option(
    '--min-clip',
    default=1,
    type=int,
    help='Require clip to extend past ref contig end by at least N bases. Default: 1',
)
@click.option(
    '--max-break',
    default=50,
    type=int,
    help='Tolerate max N unaligned bases before contig end. Default: 50',
)
@click.option(
    '--min-anchor',
    default=100,
    type=int,
    help='Minimum anchored alignment length required (default: 100).',
)
@click.option(
    '--min-mapq',
    default=0,
    type=int,
    help='Minimum mapping quality required (default: 0).',
)
@click.option(
    '--keep-secondary',
    is_flag=True,
    help='If set, include secondary alignments in output. Default: Off (exclude secondary alignments).',
)
@click.option(
    '--include-stats',
    is_flag=True,
    help='Include mapping quality, clip length, and motif counts in FASTA headers.',
)
@click.option(
    '--count-motifs',
    type=str,
    help='Comma-delimited motif sequences to count in overhang regions (e.g., "TTAGGG,CCCTAA").',
)
@click.option(
    '--fuzzy-count',
    is_flag=True,
    help='Use fuzzy motif matching allowing ±1 character variation when counting motifs.',
)
@click.option(
    '--buffer-size',
    default=1000,
    type=int,
    help='Number of sequences to buffer before writing (default: 1000).',
)
@click.option(
    '--output-format',
    default='fasta',
    type=click.Choice(['fasta', 'fastq']),
    help='Output format for extracted sequences (default: fasta).',
)
@click.option(
    '--report-stats',
    is_flag=True,
    help='Write extraction statistics to file in output directory.',
)
@click.option(
    '--no-mask-overhangs',
    is_flag=True,
    help='Do not convert overhang sequences to lowercase.',
)
@click.option(
    '--log-level',
    default='INFO',
    type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
    help='Logging level (default: INFO).',
)
@click.pass_context
def extract_cmd(
    ctx,
    samfile,
    ref_idx,
    prefix,
    extract_dir,
    min_clip,
    max_break,
    min_anchor,
    min_mapq,
    keep_secondary,
    include_stats,
    count_motifs,
    fuzzy_count,
    buffer_size,
    output_format,
    report_stats,
    no_mask_overhangs,
    log_level,
):
    """
    Extract overhanging reads for each end of each reference contig.

    Read SAM alignments and extracts soft-clipped sequences that
    extend beyond contig ends, with advanced filtering, motif analysis, and
    comprehensive statistics reporting. Overhang sequences are always written
    to output files organized by contig and end type.

    Parameters
    ----------
    ctx : click.Context
        Click context object.
    samfile : str
        Path to SAM/BAM file or '-' for stdin.
    ref_idx : str
        Path to reference index file (.fai).
    prefix : str
        Prefix for output filenames.
    extract_dir : str
        Directory for output files.
    min_clip : int
        Minimum clip length required.
    max_break : int
        Maximum gap from contig end to allow.
    min_anchor : int
        Minimum anchored alignment length required.
    min_mapq : int
        Minimum mapping quality required.
    keep_secondary : bool
        If True, include secondary alignments in output.
    include_stats : bool
        Include statistics in FASTA headers.
    count_motifs : str
        Comma-delimited motif sequences to count.
    fuzzy_count : bool
        Use fuzzy motif matching.
    buffer_size : int
        I/O buffer size for writing.
    output_format : str
        Output format ('fasta' or 'fastq').
    report_stats : bool
        Write extraction statistics to file in output directory.
    no_mask_overhangs : bool
        Disable overhang sequence masking.
    log_level : str
        Logging verbosity level.

    Examples
    --------

    # Basic extraction to current directory
    teloclip extract --ref-idx ref.fa.fai input.sam

    # Extract with motif analysis and statistics
    teloclip extract --ref-idx ref.fa.fai --include-stats \\
        --count-motifs TTAGGG,CCCTAA --report-stats input.sam

    # Extract with quality filtering and custom output
    teloclip extract --ref-idx ref.fa.fai \\
        --extract-dir overhangs/ --prefix sample1 --min-mapq 20 \\
        --min-anchor 1000 --output-format fastq input.sam

    # Read from stdin with fuzzy motif matching
    samtools view -h input.bam | teloclip extract --ref-idx ref.fa.fai \\
        --count-motifs TTAGGG --fuzzy-count
    """

    # Initialize logging for this command
    init_logging(log_level)

    try:
        # Load reference contig info
        logging.info(f'Loading reference contig info from: {ref_idx}')
        contig_info = read_fai(ref_idx)
        logging.info(f'Loaded {len(contig_info)} contigs from reference')

        # Prepare motif patterns if specified
        motif_patterns: Dict[str, str] = {}
        if count_motifs:
            # Parse comma-delimited motifs
            raw_motifs = [
                motif.strip().upper()
                for motif in count_motifs.split(',')
                if motif.strip()
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
                    if fuzzy_count:
                        pattern = make_fuzzy_motif_regex(motif)
                        pattern_name = f'{motif}_fuzzy'
                        logging.debug(f'Created fuzzy pattern for {motif}: {pattern}')
                    else:
                        pattern = make_motif_regex(motif)
                        pattern_name = motif
                        logging.debug(f'Created exact pattern for {motif}: {pattern}')
                    motif_patterns[pattern_name] = pattern

                logging.info(
                    f'Created {len(motif_patterns)} motif patterns for analysis'
                )
                if fuzzy_count:
                    logging.info(
                        'Using fuzzy matching (±1 character variation) for motif counting'
                    )

        # Initialize statistics tracker
        stats = ExtractionStats()

        # Track whether to exclude secondary alignments
        exclude_secondary = not keep_secondary

        # Create enhanced streaming filter
        logging.info('Processing alignments. Searching for overhangs.')
        alignments = EnhancedStreamingSamFilter(
            samfile=samfile,
            contigs=contig_info,
            max_break=max_break,
            min_clip=min_clip,
            min_anchor=min_anchor,
            min_mapq=min_mapq,
            motif_patterns=motif_patterns,
            stats=stats,
            exclude_secondary=exclude_secondary,
        )

        logging.info('Writing overhang reads by contig.')

        # Do not mask overhangs if output type is FASTQ
        if output_format == 'fastq':
            no_mask_overhangs = True
            logging.info('Disabling overhang masking for FASTQ output.')

        # Use enhanced extraction function
        final_stats = enhanced_streaming_split_by_contig(
            alignments=alignments,
            output_dir=extract_dir,
            prefix=prefix,
            output_format=output_format,
            buffer_size=buffer_size,
            include_stats=include_stats,
            mask_overhangs=not no_mask_overhangs,
            existing_stats=stats,
            use_sam_attributes=(include_stats and output_format == 'fastq'),
        )

        # Log comprehensive exclusion summary
        final_stats.log_exclusion_summary()

        # Generate and output statistics report
        reference_contigs = set(contig_info.keys())
        if report_stats:
            report_content = final_stats.generate_report(reference_contigs)

            # Build stats filename using extract_dir and prefix
            stats_dir = extract_dir if extract_dir else '.'
            stats_filename = f'{prefix}_stats.txt' if prefix else 'teloclip_stats.txt'
            stats_path = Path(stats_dir) / stats_filename

            logging.info(f'Writing statistics report to {stats_path}')
            # Ensure directory exists
            stats_path.parent.mkdir(parents=True, exist_ok=True)
            with open(stats_path, 'w') as f:
                f.write(report_content)

        # Final summary
        total_overhangs = final_stats.left_overhangs + final_stats.right_overhangs
        contigs_with_zero_overhangs = len(reference_contigs) - len(
            final_stats.contigs_with_overhangs
        )
        logging.info(
            f'Extraction complete: processed {final_stats.total_sam_lines} SAM lines, '
            f'found {total_overhangs} overhangs from {len(final_stats.contigs_with_overhangs)} contigs'
        )
        logging.info(
            f'Contig summary: {len(reference_contigs)} total, '
            f'{len(final_stats.contigs_with_overhangs)} with overhangs, '
            f'{contigs_with_zero_overhangs} with no overhangs'
        )

        # Show filter statistics if any alignments were filtered
        if final_stats.total_filtered > 0:
            logging.info('Filtering summary:')
            for filter_type, count in final_stats.filter_counts.items():
                if count > 0:
                    percentage = (
                        (count / final_stats.total_sam_lines) * 100
                        if final_stats.total_sam_lines > 0
                        else 0
                    )
                    logging.info(
                        f'  - {filter_type.replace("_", " ").title()}: {count} ({percentage:.1f}%)'
                    )

        if motif_patterns and final_stats.motif_matches:
            total_motif_matches = sum(final_stats.motif_matches.values())
            logging.info(
                f'Motif analysis: found {total_motif_matches} motif matches across all patterns'
            )

    except FileNotFoundError as e:
        logging.error(f'File not found: {e}')
        raise click.ClickException(str(e)) from e
    except Exception as e:
        logging.error(f'Error during extract operation: {e}')
        raise click.ClickException(str(e)) from e
