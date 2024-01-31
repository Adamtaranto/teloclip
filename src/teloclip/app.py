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

from teloclip._version import __version__
from teloclip.samops import processSamlines
from teloclip.seqops import read_fai, addRevComplement, crunchHomopolymers
from teloclip.logs import init_logging

import argparse
import logging
import sys


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
        help="Tolerate max N unaligned bases at contig ends.",
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
        help="If set tolerate +/- 1 variation in motif homopolymer runs \
            i.e. TTAGGG -> T{1,3}AG{2,4}. Default: Off",
    )
    parser.add_argument(
        "--noPoly",
        default=False,
        action="store_true",
        help='If set collapse homopolymer tracks within motifs before searching overhangs. \
                        i.e. "TTAGGGTTAGGGTTAGGGTTAGGGTTAGGG" -> "TAGTAGTAGTAGTAG". \
                        Useful for PacBio or ONP long reads homopolymer length errors. Default: Off',
    )
    parser.add_argument(
        "--matchAny",
        default=False,
        action="store_true",
        help="If set motif match may occur in unclipped region of alignment.",
    )
    parser.add_argument(
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
        if args.noPoly:
            motifList = crunchHomopolymers(motifList)
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
        matchAnywhere=args.matchAny,
        maxBreak=args.maxBreak,
        minClip=args.minClip,
        noRev=args.noRev,
        fuzzy=args.fuzzy,
    )
