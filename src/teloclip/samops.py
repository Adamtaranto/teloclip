import logging
import re
import sys

from teloclip.seqops import isMotifInClip


def processSamlines(
    samfile,
    ContigDict,
    motifList=[],
    matchAnywhere=False,
    maxBreak=0,
    minClip=1,
    noRev=False,
    fuzzy=False,
    minRepeats=1,
):
    # SAM line index keys
    SAM_QNAME = 0
    SAM_RNAME = 2
    SAM_POS = 3
    SAM_CIGAR = 5
    SAM_SEQ = 9

    # Start counters
    bothCount = 0
    keepCount = 0
    motifCount = 0
    removeCount = 0
    samlineCount = 0

    # Read SAM from stdin
    for line in samfile:
        keepLine = False
        leftClip = False
        rightClip = False
        # Write headers to stdout
        if line[0][0] == "@":
            sys.stdout.write(line)
            continue
        samlineCount += 1
        samline = line.split("\t")
        # Check if line contains soft-clip and no hard-clipping.
        if "S" in samline[SAM_CIGAR] and not "H" in samline[SAM_CIGAR]:
            # Get length of left and right overhangs
            leftClipLen, rightClipLen = checkClips(samline[SAM_CIGAR])
            alnLen = lenCIGAR(samline[SAM_CIGAR])
            # Check for left overhang
            if leftClipLen:
                if (int(samline[SAM_POS]) <= maxBreak) and (
                    leftClipLen >= (int(samline[SAM_POS]) + minClip)
                ):
                    # Overhang is on contig left
                    keepLine = True
                    leftClip = True
                    keepCount += 1
            # Check for right overhang
            if rightClipLen:
                try:
                    ContigLen = ContigDict[str(samline[SAM_RNAME])]
                except:
                    sys.exit(
                        "Reference sequence not found in FAI file: "
                        + str(samline[SAM_RNAME])
                    )
                alnEnd = int(samline[SAM_POS]) + alnLen
                # Check if overhang is on contig right end
                if ((ContigLen - alnEnd) <= maxBreak) and (
                    alnEnd + rightClipLen >= ContigLen + 1
                ):
                    rightClip = True
                    # Check if already found left OH
                    if not keepLine:
                        keepLine = True
                        keepCount += 1
                    else:
                        # Print to stderr
                        logging.info(
                            str(samline[SAM_QNAME])
                            + " overhang on both ends of "
                            + str(samline[SAM_RNAME])
                        )
                        bothCount += 1
            # Optional check for Telomeric repeat motifs
            if motifList and keepLine and matchAnywhere:
                # if noPoly:
                #    if any(
                #        s in crunchHomopolymers([samline[SAM_SEQ]])[0]
                #        for s in motifList
                #    ):
                #        sys.stdout.write(line)
                #        motifCount += 1
                #    else:
                #        removeCount += 1
                # else:
                if any(
                    s in samline[SAM_SEQ] for s in motifList
                ):  # TODO: Mod to allow regex pattern search
                    sys.stdout.write(line)
                    motifCount += 1
                else:
                    removeCount += 1
            elif motifList and keepLine:
                if isMotifInClip(
                    samline,
                    motifList,
                    leftClip,
                    rightClip,
                    leftClipLen,
                    rightClipLen,
                ):
                    sys.stdout.write(line)
                    motifCount += 1
                else:
                    removeCount += 1
            elif keepLine:
                sys.stdout.write(line)
            else:
                removeCount += 1
    if motifList:
        logging.info(
            f"Processed {samlineCount} SAM records.\n"
            f"Found {keepCount} alignments soft-clipped at contig ends.\n"
            f"Output {motifCount} alignments containing motif matches.\n"
            f"Discarded {removeCount} terminal alignments after filtering."
        )
    else:
        logging.info(
            f"Processed {samlineCount} SAM records.\n"
            f"Found {keepCount} alignments soft-clipped at contig ends.\n"
            f"Found {bothCount} alignments spanning entire contigs.\n"
            f"Discarded {removeCount} terminal alignments after filtering."
        )


def splitCIGAR(SAM_CIGAR):
    """
    Split CIGAR string into list of tuples with format (len,operator)
    """
    CIGARlist = list()
    for x in re.findall("[0-9]*[A-Z|=]", SAM_CIGAR):
        CIGARlist.append((int(re.findall("[0-9]*", x)[0]), re.findall("[A-Z]|=", x)[0]))
    # 174M76S --> [(174,M),(76,S)]
    # 96S154M --> [(96,S),(154,M)]
    return CIGARlist


def checkClips(SAM_CIGAR):
    """
    Get lengths of soft-clipped blocks from either end of an alignment given a CIGAR string.
    """
    leftClipLen = None
    rightClipLen = None
    CIGARlist = splitCIGAR(SAM_CIGAR)
    # Check if first segment is soft-clipped
    if CIGARlist[0][1] == "S":
        leftClipLen = int(CIGARlist[0][0])
    # Check if last segment is soft-clipped
    if CIGARlist[-1][1] == "S":
        rightClipLen = int(CIGARlist[-1][0])
    return (leftClipLen, rightClipLen)


def lenCIGAR(SAM_CIGAR):
    """
    Calculate length of alignment in reference sequence as sum of
    match, read-deletion, splice, mismatch, and read-match block values.
    Ignore read-insertions, padding, hard and soft clip blocks.
    """
    alnLen = 0
    CIGARlist = splitCIGAR(SAM_CIGAR)
    for x in CIGARlist:  # i.e. = [(174,M),(76,S)]
        if x[1] in set(["D", "M", "N", "X", "="]):
            alnLen += x[0]
    # Ignore operators in set('P','H','S','I')
    return alnLen


def StreamingSamFilter(samfile=None, contigs=None, maxBreak=50, minClip=1):
    """Rewrite loadSam() as generator."""
    SAM_QNAME = 0
    SAM_RNAME = 2
    SAM_POS = 3
    SAM_CIGAR = 5
    SAM_SEQ = 9
    # Read sam from stdin
    for line in samfile:
        # Skip header rows
        if line[0][0] == "@":
            continue
        samline = line.split("\t")
        # Check that aln contains soft clipping
        if "S" in samline[SAM_CIGAR] and not "H" in samline[SAM_CIGAR]:
            # Get L/R clip lengths
            leftClipLen, rightClipLen = checkClips(samline[SAM_CIGAR])
            alnLen = lenCIGAR(samline[SAM_CIGAR])
            # Check for left overhang
            if leftClipLen:
                if (int(samline[SAM_POS]) <= maxBreak) and (
                    leftClipLen >= (int(samline[SAM_POS]) + minClip)
                ):
                    # Overhang is on contig left
                    alnEnd = int(samline[SAM_POS]) + alnLen
                    try:
                        yield (
                            (
                                samline[SAM_POS],
                                alnEnd,
                                leftClipLen,
                                samline[SAM_SEQ],
                                samline[SAM_QNAME],
                                samline[SAM_RNAME],
                                "L",
                            )
                        )
                    except:
                        logging.warning(
                            "Reference sequence not found in FAI file: "
                            + str(samline[SAM_RNAME])
                        )
            # Check for right overhang
            if rightClipLen:
                try:
                    ContigLen = contigs[str(samline[SAM_RNAME])]
                except:
                    logging.warning(
                        "Reference sequence not found in FAI file: "
                        + str(samline[SAM_RNAME])
                    )
                alnEnd = int(samline[SAM_POS]) + alnLen
                # Check if overhang is on contig right end
                if ((ContigLen - alnEnd) <= maxBreak) and (
                    alnEnd + rightClipLen >= ContigLen + 1
                ):
                    yield (
                        (
                            samline[SAM_POS],
                            alnEnd,
                            rightClipLen,
                            samline[SAM_SEQ],
                            samline[SAM_QNAME],
                            samline[SAM_RNAME],
                            "R",
                        )
                    )


def SAMinfo():
    """
    Print samfile spec.
    """
    print(
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
    )


def CIGARinfo():
    """
    Print CIGAR file spec.
    """
    print(
        """
    CIGAR Operators
    --------------------
    D   Deletion; the nucleotide is present in the reference but not in the read
    H   Hard Clipping; the clipped nucleotides are not present in the read.
    I   Insertion; the nucleotide is present in the read  but not in the rference.
    M   Match; can be either an alignment match or mismatch. The nucleotide is present in the reference.
    N   Skipped region; a region of nucleotides is not present in the read
    P   Padding; padded area in the read and not in the reference
    S   Soft Clipping;  the clipped nucleotides are present in the read
    X   Read Mismatch; the nucleotide is present in the reference
    =   Read Match; the nucleotide is present in the reference

    The output order in the array is “MIDNSHP=X” followed by a field for the NM tag. If the NM tag is not present, this field will always be 0.

    M   BAM_CMATCH  0
    I   BAM_CINS    1
    D   BAM_CDEL    2
    N   BAM_CREF_SKIP   3
    S   BAM_CSOFT_CLIP  4
    H   BAM_CHARD_CLIP  5
    P   BAM_CPAD    6
    =   BAM_CEQUAL  7
    X   BAM_CDIFF   8
    B   BAM_CBACK   9
    NM  NM tag  10  
    """
    )
