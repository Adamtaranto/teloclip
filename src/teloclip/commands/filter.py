"""
Filter sub-command for teloclip CLI.

Filters SAM files for clipped alignments containing unassembled telomeric repeats.
"""

import sys

import click

from teloclip.motifs import make_fuzzy_motif_regex, make_motif_regex
from teloclip.samops import processSamlines
from teloclip.seqops import addRevComplement, read_fai


@click.command('filter')
@click.argument('samfile', type=click.File('r'), default=sys.stdin)
@click.option(
    '--ref-idx',
    required=True,
    type=click.Path(exists=True),
    help='Path to fai index for reference fasta. Index fasta using `samtools faidx FASTA`',
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
    '--motifs',
    type=str,
    help='If set keep only reads containing given motif/s from comma delimited list of strings. '
    'By default also search for reverse complement of motifs. i.e. TTAGGG,TTAAGGG will also match CCCTAA,CCCTTAA',
)
@click.option(
    '--no-rev',
    is_flag=True,
    help='If set do NOT search for reverse complement of specified motifs.',
)
@click.option(
    '--fuzzy',
    is_flag=True,
    help='If set, tolerate +/- 1 variation in motif homopolymer runs '
    'i.e. TTAGGG -> T{1,3}AG{2,4}. Default: Off',
)
@click.option(
    '-r',
    '--min-repeats',
    default=1,
    type=int,
    help='Minimum number of sequential pattern matches required for a hit to be reported. Default: 1',
)
@click.option(
    '--min-anchor',
    default=500,
    type=int,
    help='Minimum number of aligned bases (anchor) required on the non-clipped portion of the read. Default: 500',
)
@click.option(
    '--match-anywhere',
    is_flag=True,
    help='If set, motif match may occur in unclipped region of reads.',
)
@click.pass_context
def filter_cmd(
    ctx,
    samfile,
    ref_idx,
    min_clip,
    max_break,
    motifs,
    no_rev,
    fuzzy,
    min_repeats,
    min_anchor,
    match_anywhere,
):
    """
    Filter SAM file for clipped alignments containing unassembled telomeric repeats.

    Reads SAM/BAM alignments and outputs only those that are soft-clipped at contig
    ends, optionally filtering for specific motifs like telomeric repeats.

    Examples:

    \\b
    # Basic filtering for terminal soft-clipped reads
    teloclip filter --ref-idx ref.fa.fai input.sam > output.sam

    \\b
    # Filter for telomeric motifs with fuzzy matching
    teloclip filter --ref-idx ref.fa.fai --motifs TTAGGG --fuzzy input.sam

    \\b
    # Read from stdin, write to stdout
    samtools view -h input.bam | teloclip filter --ref-idx ref.fa.fai
    """
    # Fetch contig lengths
    contig_dict = read_fai(ref_idx)

    # Process motifs if provided
    motif_list = []
    if motifs:
        # Create list of motifs from comma delimited string
        motif_list = motifs.split(',')

        # Add rev comp motifs to list, make unique set.
        if not no_rev:
            motif_list = addRevComplement(motif_list)

        # Make unique set.
        motif_list = list(set(motif_list))

        # Create regex patterns for each motif
        if fuzzy:
            # Create fuzzy regex patterns for each motif
            click.echo(
                'Using option --fuzzy to find inexact motif matches. '
                'Tolerate +/- 1 base variance in motif homopolymers.',
                err=True,
            )
            motif_list = [make_fuzzy_motif_regex(motif) for motif in motif_list]
        else:
            # Create exact regex patterns for each motif
            motif_list = [make_motif_regex(motif) for motif in motif_list]

    # Process SAM lines
    processSamlines(
        samfile,
        contig_dict,
        motif_list,
        match_anywhere=match_anywhere,
        max_break=max_break,
        min_clip=min_clip,
        min_repeats=min_repeats,
        min_anchor=min_anchor,
    )
