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
    Extract overhanging reads for each end of each reference contig.

    Reads SAM/BAM alignments and extracts soft-clipped sequences that extend
    beyond contig ends, writing them to separate FASTA files for each contig end.

    Examples:

    \\b
    # Extract reads to current directory
    teloclip extract --ref-idx ref.fa.fai --extract-reads input.sam

    \\b
    # Extract with custom directory and prefix
    teloclip extract --ref-idx ref.fa.fai --extract-reads \\
        --extract-dir overhangs/ --prefix sample1 input.sam

    \\b
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
