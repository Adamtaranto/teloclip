"""
Filter sub-command for teloclip CLI.

Filters SAM files for clipped alignments containing unassembled telomeric repeats.
"""

import logging
import sys

import click

from teloclip.logs import init_logging
from teloclip.motifs import make_fuzzy_motif_regex, make_motif_regex
from teloclip.samops import processSamlines
from teloclip.seqops import addRevComplement, read_fai


@click.command(
    'filter',
    help='Filter SAM file for clipped alignments containing unassembled telomeric repeats.',
)
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
    help='Require clip to extend past ref contig end by at least N bases. Default: 1',
)
@click.option(
    '--max-break',
    default=50,
    type=int,
    help='Tolerate max N unaligned bases before contig end. Default: 50',
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
    '--keep-secondary',
    is_flag=True,
    help='If set, include secondary alignments in output. Default: Off (exclude secondary alignments).',
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
    default=100,
    type=int,
    help='Minimum number of aligned bases (anchor) required on the non-clipped portion of the read. Default: 100',
)
@click.option(
    '--match-anywhere',
    is_flag=True,
    help='If set, motif match may occur in unclipped region of reads.',
)
@click.option(
    '--log-level',
    default='INFO',
    type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR'], case_sensitive=False),
    help='Logging level (default: INFO).',
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
    keep_secondary,
    fuzzy,
    min_repeats,
    min_anchor,
    match_anywhere,
    log_level,
):
    """
    Filter SAM file for clipped alignments containing unassembled telomeric repeats.

    Reads SAM/BAM alignments and outputs only those that are soft-clipped at contig
    ends, optionally filtering for specific motifs like telomeric repeats.

    Parameters
    ----------
    ctx : click.Context
        Click context object.
    samfile : file-like
        SAM format file or stream to process.
    ref_idx : str
        Path to reference FASTA index (.fai) file.
    min_clip : int
        Minimum soft-clipping length required.
    max_break : int
        Maximum distance from contig end to consider as terminal.
    motifs : str
        Comma-separated list of DNA motifs to search for.
    no_rev : bool
        If True, do not search for reverse complement motifs.
    keep_secondary : bool
        If True, include secondary alignments in output.
    fuzzy : bool
        If True, use fuzzy matching allowing +/- 1 base variance in homopolymers.
    min_repeats : int
        Minimum number of motif repeats required for a match.
    min_anchor : int
        Minimum number of aligned bases required on non-clipped portion.
    match_anywhere : bool
        If True, allow motif matches anywhere in read, not just clipped regions.
    log_level : str
        Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).

    Examples
    --------
    Basic filtering for terminal soft-clipped reads:

    .. code-block:: bash

        teloclip filter --ref-idx ref.fa.fai input.sam > output.sam

    Filter for telomeric motifs with fuzzy matching:

    .. code-block:: bash

        teloclip filter --ref-idx ref.fa.fai --motifs TTAGGG --fuzzy input.sam

    Read from stdin, write to stdout:

    .. code-block:: bash

        samtools view -h input.bam | teloclip filter --ref-idx ref.fa.fai
    """
    # Initialize logging
    init_logging(log_level)

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

        # Log motifs being searched for
        logging.info(f'Searching for motifs: {", ".join(motif_list)}')

        # Create regex patterns for each motif
        if fuzzy:
            # Create fuzzy regex patterns for each motif
            logging.info('Tolerate +/- 1 base variance in motif homopolymers.')
            motif_list = [make_fuzzy_motif_regex(motif) for motif in motif_list]
        else:
            # Create exact regex patterns for each motif
            motif_list = [make_motif_regex(motif) for motif in motif_list]

    exclude_secondary = not keep_secondary

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
        exclude_secondary=exclude_secondary,
    )
