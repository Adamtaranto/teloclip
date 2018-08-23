#!/usr/bin/env python

from __future__ import print_function
from teloclip import __version__
import teloclip
import argparse
import sys

def log(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def mainArgs():
    parser = argparse.ArgumentParser(description='Filter SAM file for clipped alignments containing unassembled telomeric repeats.',prog='teloclip')
    # Input options
    parser.add_argument('samfile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
    parser.add_argument('--refIdx',type=str,required=True,help='Path to fai index for reference fasta. Index fasta using `samtools faidx FASTA`')
    parser.add_argument('--minClip',type=int,default=1,help='Require clip to extend past ref contig end by at least N bases.')
    parser.add_argument('--maxBreak',type=int,default=50,help='Tolerate max N unaligned bases at contig ends.')
    parser.add_argument('--motifs',type=str,default=None, help='If set keep only reads containing given motif/s from comma delimited list of strings. \
    By default also search for reverse complement of motifs. i.e. TTAGGG,TTAAGGG will also match CCCTAA,CCCTTAA')
    parser.add_argument('--norev',default=False,action='store_true', help='If set do NOT search for reverse complement of specified motifs.')
    parser.add_argument('--matchAny',default=False,action='store_true', help='If set motif match may occur in unclipped region of alignment.')
    parser.add_argument('--version', action='version',version='%(prog)s {version}'.format(version=__version__))
    args = parser.parse_args()
    return args

def main():
    # Get cmd line args
    args = mainArgs()
    # SAM line index keys
    SAM_QNAME   = 0
    SAM_RNAME   = 2
    SAM_POS     = 3
    SAM_CIGAR   = 5
    SAM_SEQ     = 9
    # Start counters
    bothCount = 0
    keepCount = 0
    motifCount = 0
    removeCount = 0
    samlineCount = 0
    # Fetch contig lens
    ContigDict = teloclip.read_fai(args.refIdx)
    # Split motifs into list
    if args.motifs:
        motifList = args.motifs.split(",")
        # Add rev comp motifs to list, make unique set.
        if not args.norev:
            motifList = teloclip.addRevComplement(motifList)
    # Read SAM from stdin
    for line in args.samfile:
        keepLine = False
        leftClip = False
        rightClip = False
        if line[0][0] == "@":
            sys.stdout.write(line)
            continue
        samlineCount += 1
        samline = line.split('\t')
        # Check if line contains soft-clip and no hard-clipping.
        if "S" in samline[SAM_CIGAR] and not "H" in samline[SAM_CIGAR]:
            leftClipLen,rightClipLen = teloclip.checkClips(samline[SAM_CIGAR])
            alnLen = teloclip.lenCIGAR(samline[SAM_CIGAR])
            # Check for left overhang
            if leftClipLen:
                if (int(samline[SAM_POS]) <= args.maxBreak) and (leftClipLen >= (int(samline[SAM_POS]) + args.minClip)):
                    # Overhang is on contig left
                    keepLine = True
                    leftClip = True
                    keepCount += 1
            # Check for right overhang
            if rightClipLen:
                try: ContigLen = ContigDict[str(samline[SAM_RNAME])]
                except: sys.exit("Reference sequence not found in FAI file: " + str(samline[SAM_RNAME]))
                alnEnd = int(samline[SAM_POS]) + alnLen 
                # Check if overhang is on contig right end
                if ((ContigLen - alnEnd) <= args.maxBreak) and (alnEnd + rightClipLen >= ContigLen +1 ):
                    rightClip = True
                    # Check if already found left OH
                    if not keepLine:
                        keepLine = True
                        keepCount += 1
                    else: 
                        # Print to stderr
                        log(str(samline[SAM_QNAME]) + " overhang on both ends of " + str(samline[SAM_RNAME]))
                        bothCount += 1
            # Optional check for Telomeric repeat motifs
            if args.motifs and keepLine and args.matchAny:
                if any(s in samline[SAM_SEQ] for s in motifList):
                    sys.stdout.write(line)
                    motifCount += 1
                else:
                    removeCount += 1
            elif args.motifs and keepLine:
                if teloclip.isClipMotif(samline,motifList,leftClip,rightClip,leftClipLen,rightClipLen):
                    sys.stdout.write(line)
                    motifCount += 1
                else:
                    removeCount += 1
            elif keepLine:  
                sys.stdout.write(line)
            else:
                removeCount += 1
    if args.motifs:
        log("Processes %s SAM records. \nFound %s alignments soft-clipped at contig ends. \nOutput %s alignments containing motif matches. \nDiscarded %s terminal alignments after filtering." % (samlineCount,keepCount,motifCount,removeCount))
    else:
        log("Processes %s SAM records. \nFound %s alignments soft-clipped at contig ends. \nFound %s alignments spanning entire contigs. \nDiscarded %s terminal alignments after filtering." % (samlineCount,keepCount,bothCount,removeCount))

"""
# SAM format
# ----------------------------------------------------------------
# 0   QNAME   String  # Query template NAME
# 1   FLAG    Int     # Bitwise FLAG
# 2   RNAME   String  # References sequence NAME
# 3   POS     Int     # 1-based leftmost mapping POSition
# 4   MAPQ    Int     # MAPping Quality
# 5   CIGAR   String  # CIGAR String
# 6   RNEXT   String  # Ref. name of the mate/next read
# 7   PNEXT   Int     # Position of the mate/next read
# 8   TLEN    Int     # Observed Template LENgth
# 9   SEQ     String  # Segment SEQuence
# 10  QUAL    String  # ASCII of Phred-scaled base QUALity+33
"""