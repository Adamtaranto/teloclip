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
from teloclip.samops import processSamlines
from teloclip.seqops import read_fai, addRevComplement


def mainArgs():
    parser = argparse.ArgumentParser(
        description="Filter SAM file for clipped alignments containing unassembled telomeric repeats.",
        prog="teloclip",
    )
    # Input options
    parser.add_argument(
        "samfile", nargs="?", type=argparse.FileType("r"), default=sys.stdin
    )
    parser.add_argument(
        "--refIdx",
        type=str,
        required=True,
        help="Path to fai index for reference fasta. Index fasta using `samtools faidx FASTA`",
    )
    parser.add_argument(
        "--minClip",
        type=int,
        default=1,
        help="Require clip to extend past ref contig end by at least N bases.",
    )
    parser.add_argument(
        "--maxBreak",
        type=int,
        default=50,
        help="Tolerate max N unaligned bases before contig end.",
    )
    parser.add_argument(
        "--motifs",
        type=str,
        default=None,
        help="If set keep only reads containing given motif/s from comma delimited list of strings. \
                        By default also search for reverse complement of motifs. i.e. TTAGGG,TTAAGGG will also match CCCTAA,CCCTTAA",
    )
    parser.add_argument(
        "--noRev",
        default=False,
        action="store_true",
        help="If set do NOT search for reverse complement of specified motifs.",
    )
    parser.add_argument(
        "--fuzzy",
        default=False,
        action="store_true",
        help="If set, tolerate +/- 1 variation in motif homopolymer runs \
            i.e. TTAGGG -> T{1,3}AG{2,4}. Default: Off",
    )
    parser.add_argument(
        "-r",
        "--min_repeats",
        default=3,
        type=int,
        help="Minimum number of sequential pattern matches required for a hit to be reported. Default: 3",
    )
    parser.add_argument(
        "--noPoly",
        default=False,
        action="store_true",
        help="WARNING: OPTION DEPRECIATED, USE --FUZZY. Default: False",
    )
    parser.add_argument(
        "--matchAnywhere",
        default=False,
        action="store_true",
        help="If set, motif match may occur in unclipped region of alignment.",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )
    args = parser.parse_args()
    return args


def main():
    # Set up logging
    init_logging()

    # Get cmd line args
    args = mainArgs()

    # Fetch contig lens
    ContigDict = read_fai(args.refIdx)

    # Split motifs into list
    if args.motifs:
        motifList = args.motifs.split(",")
        # Add rev comp motifs to list, make unique set.
        if not args.noRev:
            motifList = addRevComplement(motifList)
        if args.fuzzy:
            motifList = motifList  # TODO: Implement convert to regex
            # convert to regex pattern(motifList)
    else:
        motifList = []

    # Log depreciation warning if using --noPoly option
    if args.noPoly:
        logging.warning(
            "WARNING: Option --noPoly is depreciated and should not be used. \
            Switching to option --fuzzy instead to find inexact motif matches."
        )
        args.fuzzy = True

    if args.fuzzy:
        logging.info(
            "INFO: Using option --fuzzy to find inexact motif matches. \
                Tolerate +/- 1 base variance in motif homopolymers."
        )

    processSamlines(
        args.samfile,
        ContigDict,
        motifList,
        matchAnywhere=args.matchAnywhere,
        maxBreak=args.maxBreak,
        minClip=args.minClip,
        noRev=args.noRev,
        fuzzy=args.fuzzy,
        minRepeats=args.min_repeats,
    )
