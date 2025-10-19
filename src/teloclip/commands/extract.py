"""
Extract sub-command for teloclip CLI.

Extract overhanging reads for each end of each reference contig and write to FASTA files.
"""

import logging
import sys

import click

from teloclip.samops import StreamingSamFilter, StreamingSplitByContig
from teloclip.seqops import read_fai


@click.command('extract')
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
@click.pass_context
def extract_cmd(
    ctx, samfile, ref_idx, prefix, extract_reads, extract_dir, min_clip, max_break
):
    """
    Extract overhanging (soft-clipped) read sequences that extend beyond contig ends.

    Reads SAM/BAM alignments, identifies reads with soft-clipped segments that
    extend past either end of reference contigs (as defined in a FASTA index),
    and optionally writes those overhanging sequences to per-contig FASTA files.
    If writing is disabled, the function still consumes and processes the
    alignment stream to trigger any side effects of the underlying streaming
    pipeline.

    Parameters
    ----------
    ctx : click.Context
        Click command context passed by the CLI entry point. Used for CLI-related
        state and logging but not required for core processing.
    samfile : str or file-like
        Path to a SAM/BAM file or a file-like object/stream. Use '-' or pipe data
        in via stdin to stream SAM records.
    ref_idx : str
        Path to the reference FASTA index file (FAI). Used to obtain contig names
        and lengths to determine contig ends.
    prefix : str
        Filename prefix to use when writing output FASTA files for extracted
        overhang reads. Combined with contig identifiers to form output names.
    extract_reads : bool
        If True, extracted overhang sequences are written to FASTA files, one
        file (or set of files) per contig end. If False, alignments are still
        processed but no files are written.
    extract_dir : str
        Destination directory for output FASTA files when extract_reads is True.
        The directory will be created if necessary (behavior depends on the
        underlying implementation).
    min_clip : int
        Minimum soft-clip length (in bases) required for a clipped segment to be
        considered an overhang. Clipped segments shorter than this are ignored.
    max_break : int
        Maximum allowed internal break/gap (in bases) within an alignment used by
        the streaming filter to decide whether a clipped segment represents a
        contiguous overhang.

    Returns
    -------
    None
        This function performs streaming I/O and side effects (logging and optional
        file output) and does not return a value.

    Notes
    -----
    - Contig information is loaded from the provided FASTA index (ref_idx) to
      determine contig boundaries.
    - Alignments are streamed and filtered; when extract_reads is True, matching
      soft-clipped sequences are grouped by contig end and written as FASTA.
      When extract_reads is False, the alignment stream is still consumed to
      perform processing (but no files are produced).
    - Underlying I/O, parsing, or streaming utilities may raise exceptions
      (e.g., FileNotFoundError for missing files, parsing errors for malformed
      SAM/BAM/FAI). Those exceptions are propagated to the caller.

    Examples
    --------
    # Extract reads to current directory
    teloclip extract --ref-idx ref.fa.fai --extract-reads input.sam

    # Extract with custom directory and prefix
    teloclip extract --ref-idx ref.fa.fai --extract-reads \\
        --extract-dir overhangs/ --prefix sample1 input.sam

    # Read from stdin
    samtools view -h input.bam | teloclip extract --ref-idx ref.fa.fai --extract-reads
    """
    # Load ref contigs lengths as dict
    logging.info('Importing reference contig info from: %s', ref_idx)
    contig_info = read_fai(ref_idx)

    # Load alignments from samfile or stdin
    logging.info('Processing alignments. Searching for overhangs.')

    alignments = StreamingSamFilter(
        samfile=samfile,
        contigs=contig_info,
        max_break=max_break,
        min_clip=min_clip,
    )

    if extract_reads:
        logging.info('Writing overhang reads by contig.')
        StreamingSplitByContig(
            alignments=alignments,
            contigs=contig_info,
            prefix=prefix,
            outdir=extract_dir,
        )
    else:
        # If extract_reads is not set, we still need to consume the generator
        # to process the alignments, but we won't write anything
        for _ in alignments:
            pass
        logging.info('Processing complete. Use --extract-reads to write output files.')
