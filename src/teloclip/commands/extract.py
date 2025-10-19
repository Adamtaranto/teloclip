"""
Extract sub-command for teloclip CLI.

This refactored version provides significant performance improvements and new features:
- Memory-efficient I/O with BioPython SeqIO integration
- Rich sequence headers with mapping quality and motif statistics
- Motif analysis integration with exact and fuzzy matching
- Comprehensive statistics reporting
- Quality filtering and validation
- Support for FASTA and FASTQ output formats
"""

import logging
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
    'extract', help='Extract overhanging reads for each end of each reference contig.'
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
    '--extract-reads',
    is_flag=True,
    help='If set, write overhang reads to fasta by contig.',
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
    help='Require clip to extend past ref contig end by at least N bases.',
)
@click.option(
    '--max-break',
    default=50,
    type=int,
    help='Tolerate max N unaligned bases before contig end.',
)
@click.option(
    '--min-anchor',
    default=500,
    type=int,
    help='Minimum anchored alignment length required (default: 500).',
)
@click.option(
    '--min-mapq',
    default=0,
    type=int,
    help='Minimum mapping quality required (default: 0).',
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
    '--stats-report',
    type=click.Path(),
    help='Write extraction statistics to file. Use "-" for stdout.',
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
    extract_reads,
    extract_dir,
    min_clip,
    max_break,
    min_anchor,
    min_mapq,
    include_stats,
    count_motifs,
    fuzzy_count,
    buffer_size,
    output_format,
    stats_report,
    no_mask_overhangs,
    log_level,
):
    """
    Extract overhanging reads for each end of each reference contig.

    This enhanced version reads SAM/BAM alignments and extracts soft-clipped sequences
    that extend beyond contig ends, with advanced filtering, motif analysis, and
    comprehensive statistics reporting.

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
    extract_reads : bool
        Extract and write overhang sequences.
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
    stats_report : str
        Path for statistics report output.
    no_mask_overhangs : bool
        Disable overhang sequence masking.
    log_level : str
        Logging verbosity level.

    Examples
    --------

    # Basic extraction to current directory
    teloclip extract --ref-idx ref.fa.fai --extract-reads input.sam

    # Extract with motif analysis and statistics
    teloclip extract --ref-idx ref.fa.fai --extract-reads --include-stats \\
        --count-motifs TTAGGG,CCCTAA --stats-report stats.txt input.sam

    # Extract with quality filtering and custom output
    teloclip extract --ref-idx ref.fa.fai --extract-reads \\
        --extract-dir overhangs/ --prefix sample1 --min-mapq 20 \\
        --min-anchor 1000 --output-format fastq input.sam

    # Read from stdin with fuzzy motif matching
    samtools view -h input.bam | teloclip extract --ref-idx ref.fa.fai \\
        --extract-reads --count-motifs TTAGGG --fuzzy-count
    """
    # Initialize logging
    init_logging(level=getattr(logging, log_level.upper()))
    logger = logging.getLogger(__name__)

    try:
        # Load reference contig info
        logger.info(f'Loading reference contig info from: {ref_idx}')
        contig_info = read_fai(ref_idx)
        logger.info(f'Loaded {len(contig_info)} contigs from reference')

        # Prepare motif patterns if specified
        motif_patterns: Dict[str, str] = {}
        if count_motifs:
            # Parse comma-delimited motifs
            raw_motifs = [
                motif.strip().upper()
                for motif in count_motifs.split(',')
                if motif.strip()
            ]
            logger.info(f'Processing motif list: {", ".join(raw_motifs)}')

            # Validate motifs - must contain only A, T, G, C
            valid_bases = {'A', 'T', 'G', 'C'}
            validated_motifs = []

            for motif in raw_motifs:
                if not motif:  # Skip empty motifs
                    continue
                if not all(base in valid_bases for base in motif):
                    invalid_bases = set(motif) - valid_bases
                    logger.warning(
                        f'Skipping invalid motif "{motif}": contains invalid bases {invalid_bases}'
                    )
                    continue
                validated_motifs.append(motif)

            if not validated_motifs:
                logger.warning('No valid motifs found after validation')
            else:
                logger.info(
                    f'Validated {len(validated_motifs)} motifs: {", ".join(validated_motifs)}'
                )

                # Add reverse complements and create unique set
                all_motifs = set()
                for motif in validated_motifs:
                    all_motifs.add(motif)
                    rev_comp = revComp(motif)
                    all_motifs.add(rev_comp)
                    logger.debug(f'Motif: {motif} -> Reverse complement: {rev_comp}')

                # Convert to sorted list for consistent ordering
                unique_motifs = sorted(all_motifs)
                logger.info(
                    f'Final motif set (including reverse complements): {", ".join(unique_motifs)}'
                )

                # Create regex patterns for each unique motif
                for motif in unique_motifs:
                    if fuzzy_count:
                        pattern = make_fuzzy_motif_regex(motif)
                        pattern_name = f'{motif} (fuzzy)'
                        logger.debug(f'Created fuzzy pattern for {motif}: {pattern}')
                    else:
                        pattern = make_motif_regex(motif)
                        pattern_name = motif
                        logger.debug(f'Created exact pattern for {motif}: {pattern}')
                    motif_patterns[pattern_name] = pattern

                logger.info(
                    f'Created {len(motif_patterns)} motif patterns for analysis'
                )
                if fuzzy_count:
                    logger.info(
                        'Using fuzzy matching (±1 character variation) for motif counting'
                    )

        # Initialize statistics tracker
        stats = ExtractionStats()

        # Create enhanced streaming filter
        logger.info('Processing alignments. Searching for overhangs.')
        alignments = EnhancedStreamingSamFilter(
            samfile=samfile,
            contigs=contig_info,
            max_break=max_break,
            min_clip=min_clip,
            min_anchor=min_anchor,
            min_mapq=min_mapq,
            motif_patterns=motif_patterns,
            stats=stats,
        )

        if extract_reads:
            logger.info('Writing overhang reads by contig.')

            # Use enhanced extraction function
            final_stats = enhanced_streaming_split_by_contig(
                alignments=alignments,
                output_dir=extract_dir,
                prefix=prefix,
                output_format=output_format,
                buffer_size=buffer_size,
                include_stats=include_stats,
                mask_overhangs=not no_mask_overhangs,
            )

        else:
            # Process alignments without writing files (for statistics only)
            for alignment in alignments:
                stats.record_alignment(
                    contig_name=alignment['contig_name'],
                    is_left=(alignment['end'] == 'L'),
                    motif_counts=alignment.get('motif_counts'),
                )
            final_stats = stats
            logger.info(
                'Processing complete. Use --extract-reads to write output files.'
            )

        # Generate and output statistics report
        if stats_report or not extract_reads:
            report_content = final_stats.generate_report()

            if stats_report:
                if stats_report == '-':
                    logger.info('Writing statistics report to stdout')
                    print(report_content)
                else:
                    logger.info(f'Writing statistics report to {stats_report}')
                    with open(stats_report, 'w') as f:
                        f.write(report_content)
            elif not extract_reads:
                # Always show stats if not extracting files
                print(report_content)

        # Final summary
        logger.info(
            f'Extraction complete: processed {final_stats.total_alignments} alignments, '
            f'found {final_stats.left_overhangs + final_stats.right_overhangs} overhangs '
            f'from {len(final_stats.contigs_processed)} contigs'
        )

        if motif_patterns and final_stats.motif_matches:
            total_motif_matches = sum(final_stats.motif_matches.values())
            logger.info(
                f'Motif analysis: found {total_motif_matches} motif matches across all patterns'
            )

    except FileNotFoundError as e:
        logger.error(f'File not found: {e}')
        raise click.ClickException(str(e)) from e
    except Exception as e:
        logger.error(f'Error during extract operation: {e}')
        raise click.ClickException(str(e)) from e
