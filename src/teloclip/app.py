"""
                                       ██████╗ ██╗      ██╗ ██████╗
                                      ██╔════╝ ██║      ██║ ██╔══██╗
                                      ██║      ██║      ██║ ██████╔╝
████████╗ ███████╗ ██╗       ██████╗  ██║      ██║      ██║ ██╔═══╝
╚══██╔══╝ ██╔════╝ ██║      ██╔═══██╗ ╚██████╗ ███████╗ ██║ ██║
   ██║    █████╗   ██║      ██║   ██║  ╚═════╝ ╚══════╝ ╚═╝ ╚═╝
   ██║    ██╔══╝   ██║      ██║   ██║
   ██║    ███████╗ ███████╗ ╚██████╔╝
   ╚═╝    ╚══════╝ ╚══════╝  ╚═════╝

A tool for the recovery of unassembled telomeres from soft-clipped read alignments.

"""

import argparse
import logging
import sys

from teloclip._version import __version__
from teloclip.logs import init_logging
from teloclip.motifs import make_fuzzy_motif_regex, make_motif_regex
from teloclip.samops import processSamlines
from teloclip.seqops import addRevComplement, read_fai


def mainArgs():
    parser = argparse.ArgumentParser(
        description='Filter SAM file for clipped alignments containing unassembled telomeric repeats.',
        prog='teloclip',
    )
    # Input options
    parser.add_argument(
        'samfile', nargs='?', type=argparse.FileType('r'), default=sys.stdin
    )
    parser.add_argument(
        '--ref-idx',
        dest='ref_idx',
        type=str,
        required=True,
        help='Path to fai index for reference fasta. Index fasta using `samtools faidx FASTA`',
    )
    parser.add_argument(
        '--min-clip',
        dest='min_clip',
        type=int,
        default=1,
        help='Require clip to extend past ref contig end by at least N bases.',
    )
    parser.add_argument(
        '--max-break',
        dest='max_break',
        type=int,
        default=50,
        help='Tolerate max N unaligned bases before contig end.',
    )
    parser.add_argument(
        '--motifs',
        type=str,
        default=None,
        help='If set keep only reads containing given motif/s from comma delimited list of strings. '
        'By default also search for reverse complement of motifs. i.e. TTAGGG,TTAAGGG will also match CCCTAA,CCCTTAA',
    )
    parser.add_argument(
        '--no-rev',
        dest='no_rev',
        default=False,
        action='store_true',
        help='If set do NOT search for reverse complement of specified motifs.',
    )
    parser.add_argument(
        '--fuzzy',
        default=False,
        action='store_true',
        help='If set, tolerate +/- 1 variation in motif homopolymer runs '
        'i.e. TTAGGG -> T{1,3}AG{2,4}. Default: Off',
    )
    parser.add_argument(
        '-r',
        '--min-repeats',
        dest='min_repeats',
        default=1,
        type=int,
        help='Minimum number of sequential pattern matches required for a hit to be reported. Default: 3',
    )
    parser.add_argument(
        '--min-anchor',
        dest='min_anchor',
        default=500,
        type=int,
        help='Minimum number of aligned bases (anchor) required on the non-clipped portion of the read. Default: 500',
    )
    parser.add_argument(
        '--match-anywhere',
        dest='match_anywhere',
        default=False,
        action='store_true',
        help='If set, motif match may occur in unclipped region of reads.',
    )
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__),
    )
    args = parser.parse_args()

    return args


def main():
    # Set up logging
    init_logging()

    # Get cmd line args
    args = mainArgs()

    # Fetch contig lens
    contig_dict = read_fai(args.ref_idx)

    if args.motifs:
        # Create list of motifs from comma delimited string
        motifList = args.motifs.split(',')

        # Add rev comp motifs to list, make unique set.
        if not args.no_rev:
            motifList = addRevComplement(motifList)

        # Make unique set.
        motifList = list(set(motifList))

        # Create regex patterns for each motif
        if args.fuzzy:
            # Create fuzzy regex patterns for each motif
            logging.info(
                'Using option --fuzzy to find inexact motif matches. Tolerate +/- 1 base variance in motif homopolymers.'
            )
            motifList = [make_fuzzy_motif_regex(motif) for motif in motifList]
        else:
            # Create exact regex patterns for each motif
            motifList = [make_motif_regex(motif) for motif in motifList]
    else:
        motifList = []

    # Process SAM lines

    processSamlines(
        args.samfile,
        contig_dict,
        motifList,
        match_anywhere=args.match_anywhere,
        max_break=args.max_break,
        min_clip=args.min_clip,
        min_repeats=args.min_repeats,
        min_anchor=args.min_anchor,
    )
