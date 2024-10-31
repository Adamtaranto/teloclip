import argparse
import logging
import os
import sys

from teloclip._version import __version__
from teloclip.logs import init_logging
from teloclip.samops import StreamingSamFilter
from teloclip.seqops import writefasta, read_fai

"""
Module under development. Currently splits overhang reads of each contig end into separate fata files. 
Users can then assemble these reads externally and realign to the ref contig to repair missing telomeres.
"""


def mainArgs():
    parser = argparse.ArgumentParser(
        description="Extract overhanging reads for each end of each reference contig. Write to fasta.",
        prog="teloclip-extract",
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
        "--prefix",
        type=str,
        default=None,
        required=False,
        help="Use this prefix for output files. Default: None.",
    )

    parser.add_argument(
        "--extractReads",
        default=False,
        action="store_true",
        help="If set, write overhang reads to fasta by contig.",
    )

    parser.add_argument(
        "--extractDir",
        type=str,
        default=None,
        required=False,
        help="Write extracted reads to this directory. Default: cwd.",
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

    # Version info
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )

    args = parser.parse_args()

    return args


def StreamingSplitByContig(alignments=None, contigs=None, prefix=None, outdir=None):
    """
    Takes alignment line summaries tagged with overhang information from StreamingSamFilter
    Writes reads into output files for each end of each contig with at least one alignment found.
    """
    # Note: This method opens and closes the output fasta files for EVERY read processed.
    # Probably inefficient but avoids reading everything into memory
    if outdir:
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
    else:
        outdir = os.getcwd()
    # Output file tracker
    outpaths = list()
    # Counter
    readCount = 0
    # For each aligned read
    for read in alignments:
        readCount += 1
        # Log update every 10K reads
        if not readCount % 10000:
            logging.info("Alignments processed: %s" % str(readCount))
        # Check if alignment is at right end overhang
        if read[6] == "R":
            # i.e [(alnStart,alnEnd,rightClipLen,readSeq,readName,contigName,'R')]
            if prefix:
                base = "_".join([prefix, str(read[5])])
            else:
                base = str(read[5])
            # Set output file path
            outfileR = os.path.join(outdir, "_".join([base, "R"]) + ".fasta")
            # Check if need to create new outfile or append to existing file.
            if outfileR not in outpaths:
                logging.info("Creating new outfile: %s" % str(outfileR))
                outpaths.append(outfileR)
                # Output reads aligned to right end of contig
                with open(outfileR, "w") as fileR:
                    # Convert overhanging segment of read to lowercase
                    masked = read[3][: -read[2]] + read[3][-read[2] :].lower()
                    # Write masked read to fasta
                    writefasta(fileR, str(read[4]), masked)
            else:
                # If file already exists
                # Output reads aligned to right end of contig
                with open(outfileR, "a") as fileR:
                    # Convert overhanging segment of read to lowercase
                    masked = read[3][: -read[2]] + read[3][-read[2] :].lower()
                    # Write masked read to fasta
                    writefasta(fileR, str(read[4]), masked)
        # Check if alignment is a left end overhang
        elif read[6] == "L":
            # i.e [(alnStart,alnEnd,leftClipLen,readSeq,readname,contigname,'L')]
            if prefix:
                base = "_".join([prefix, str(read[5])])
            else:
                base = str(read[5])
            # Set output file path
            outfileL = os.path.join(outdir, "_".join([base, "L"]) + ".fasta")
            # Check if need to create new outfile or append to existing file.
            if outfileL not in outpaths:
                logging.info("Creating new outfile: %s" % str(outfileL))
                outpaths.append(outfileL)
                # Output reads aligned to left end of contig
                with open(outfileL, "w") as fileL:
                    # Convert overhanging segment of read to lowercase
                    masked = read[3][: read[2]].lower() + read[3][read[2] :]
                    # Write masked read to fasta
                    writefasta(fileL, str(read[4]), masked)
            else:
                # If file already exists
                # Output reads aligned to left end of contig
                with open(outfileL, "a") as fileL:
                    # Convert overhanging segment of read to lowercase
                    masked = read[3][: read[2]].lower() + read[3][read[2] :]
                    # Write masked read to fasta
                    writefasta(fileL, str(read[4]), masked)

    logging.info("Total alignments processed: %s" % str(readCount))


def main():
    # Set up logging
    init_logging()

    # Get cmd line args
    args = mainArgs()

    # Load ref contigs lengths as dict
    logging.info("Importing reference contig info from: %s " % str(args.refIdx))
    contigInfo = read_fai(args.refIdx)

    # Load alignments from samfile or stdin
    # {'contig':{"L":[(alnStart,alnEnd,leftClipLen,readSeq,readname)],"R":[(alnStart,alnEnd,rightClipLen,readSeq,readname)]}}
    logging.info("Processing alignments. Searching for overhangs.")

    alignments = StreamingSamFilter(
        samfile=args.samfile,
        contigs=contigInfo,
        maxBreak=args.maxBreak,
        minClip=args.minClip,
    )

    if args.extractReads:
        logging.info("Writing overhang reads by contig.")
        # splitbycontig(alignments=alignments,contigs=contigInfo,prefix=args.prefix,outdir=args.extractDir)
        StreamingSplitByContig(
            alignments=alignments,
            contigs=contigInfo,
            prefix=args.prefix,
            outdir=args.extractDir,
        )

    sys.exit(0)
