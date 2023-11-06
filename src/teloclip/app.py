#!/usr/bin/env python

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

from __future__ import print_function

from teloclip._version import __version__
from teloclip.samops import processSamlines
from teloclip.seqops import read_fai, addRevComplement, crunchHomopolymers

import argparse
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

    processSamlines(
        args.samfile,
        ContigDict,
        motifList,
        matchAnywhere=args.matchAny,
        maxBreak=args.maxBreak,
        minClip=args.minClip,
        noPoly=args.noPoly,
        noRev=args.noRev,
    )
