import logging
import re
import sys
import os

from teloclip.motifs import check_sequence_for_patterns
from teloclip.seqops import isMotifInClip, writefasta

# TODO:
"""
- Create a samline class to hold sam line data
- Create a CIGAR class based on https://github.com/brentp/cigar/blob/master/cigar.py
- Sam fields as attributes of samline class
- Calculate clip regions as attributes of samline class
- Methods to check if left or right clip is terminal
- Method to check if motif in clips
- Method to print samline
- Note if left or right terminal softclip regions present
- Add and use min_anchor threshold (default 500bp) to ensure that the alignment is anchored to the reference
- Split left and right clip checks into helper functions
- Add function: Count Cigar len to left / right of softclip 20S20M2D10M2S = 30 bases of read aligned to ref
    - Check that aligned length on ref is >= min_anchor
"""


def processSamlines(
    samfile,
    contig_dict,
    motif_list=None,
    match_anywhere=False,
    max_break=0,
    min_clip=1,
    min_repeats=1,
    min_anchor=500,
):
    # SAM line index keys
    if motif_list is None:
        motif_list = []
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
    anchorFilteredCount = 0

    # Read SAM from stdin
    for line in samfile:
        keepLine = False
        leftClip = False
        rightClip = False
        # Write headers to stdout
        if line[0][0] == '@':
            sys.stdout.write(line)
            continue
        samlineCount += 1
        samline = line.split('\t')
        # Check if line contains soft-clip and no hard-clipping.
        if 'S' in samline[SAM_CIGAR] and 'H' not in samline[SAM_CIGAR]:
            # Check if alignment meets minimum anchor requirement
            if not validate_min_anchor(samline[SAM_CIGAR], min_anchor):
                anchorFilteredCount += 1
                removeCount += 1
                continue

            # Get length of left and right overhangs
            leftClipLen, rightClipLen = checkClips(samline[SAM_CIGAR])
            alnLen = lenCIGAR(samline[SAM_CIGAR])
            # Check for left overhang
            if leftClipLen:
                if (int(samline[SAM_POS]) <= max_break) and (
                    leftClipLen >= (int(samline[SAM_POS]) + min_clip)
                ):
                    # Overhang is on contig left
                    keepLine = True
                    leftClip = True
                    keepCount += 1
            # Check for right overhang
            if rightClipLen:
                try:
                    ContigLen = contig_dict[str(samline[SAM_RNAME])]
                except KeyError:
                    sys.exit(
                        'Reference sequence not found in FAI file: '
                        + str(samline[SAM_RNAME])
                    )
                alnEnd = int(samline[SAM_POS]) + alnLen
                # Check if overhang is on contig right end
                if ((ContigLen - alnEnd) <= max_break) and (
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
                            + ' overhang on both ends of '
                            + str(samline[SAM_RNAME])
                        )
                        bothCount += 1
            # Optional check for Telomeric repeat motifs
            if motif_list and keepLine and match_anywhere:
                if check_sequence_for_patterns(
                    samline[SAM_SEQ], motif_list, min_repeats
                ):
                    sys.stdout.write(line)
                    motifCount += 1
                else:
                    removeCount += 1
            elif motif_list and keepLine:
                if isMotifInClip(
                    samline,
                    motif_list,
                    leftClip,
                    rightClip,
                    leftClipLen,
                    rightClipLen,
                    min_repeats,
                ):
                    sys.stdout.write(line)
                    motifCount += 1
                else:
                    removeCount += 1
            elif keepLine:
                sys.stdout.write(line)
            else:
                removeCount += 1
    if motif_list:
        logging.info(
            f'Processed {samlineCount} SAM records.\n'
            f'Found {keepCount} alignments soft-clipped at contig ends.\n'
            f'Filtered {anchorFilteredCount} alignments below min_anchor threshold ({min_anchor}bp).\n'
            f'Output {motifCount} alignments containing motif matches.\n'
            f'Discarded {removeCount} total alignments after all filtering.'
        )
    else:
        logging.info(
            f'Processed {samlineCount} SAM records.\n'
            f'Found {keepCount} alignments soft-clipped at contig ends.\n'
            f'Filtered {anchorFilteredCount} alignments below min_anchor threshold ({min_anchor}bp).\n'
            f'Found {bothCount} alignments spanning entire contigs.\n'
            f'Discarded {removeCount} total alignments after all filtering.'
        )


def splitCIGAR(SAM_CIGAR):
    """
    Split CIGAR string into list of tuples with format (len,operator)
    """
    CIGARlist = []
    for x in re.findall('[0-9]*[A-Z|=]', SAM_CIGAR):
        CIGARlist.append((int(re.findall('[0-9]*', x)[0]), re.findall('[A-Z]|=', x)[0]))
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
    if CIGARlist[0][1] == 'S':
        leftClipLen = int(CIGARlist[0][0])
    # Check if last segment is soft-clipped
    if CIGARlist[-1][1] == 'S':
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
        if x[1] in {'D', 'M', 'N', 'X', '='}:
            alnLen += x[0]
    # Ignore operators in set('P','H','S','I')
    return alnLen


def calculate_aligned_bases(cigar_string):
    """
    Calculate the number of bases that are actually aligned/matched to the reference.
    Only counts M (match), = (sequence match), and X (mismatch) operations.
    Excludes deletions (D), insertions (I), soft clips (S), hard clips (H),
    padding (P), and splicing (N) operations.

    This is used for min_anchor validation to ensure sufficient anchoring
    of the read to the reference sequence.

    Args:
        cigar_string (str): CIGAR string from SAM alignment

    Returns:
        int: Number of aligned bases (M + = + X operations only)
    """
    aligned_bases = 0
    cigar_list = splitCIGAR(cigar_string)
    for length, operation in cigar_list:
        if operation in {'M', '=', 'X'}:
            aligned_bases += length
    return aligned_bases


def validate_min_anchor(cigar_string, min_anchor):
    """
    Validate that an alignment has sufficient anchored bases to meet min_anchor requirement.

    Args:
        cigar_string (str): CIGAR string from SAM alignment
        min_anchor (int): Minimum number of aligned bases required

    Returns:
        bool: True if alignment meets min_anchor requirement, False otherwise
    """
    aligned_bases = calculate_aligned_bases(cigar_string)
    return aligned_bases >= min_anchor


def StreamingSamFilter(samfile=None, contigs=None, max_break=50, min_clip=1):
    """Rewrite loadSam() as generator."""
    SAM_QNAME = 0
    SAM_RNAME = 2
    SAM_POS = 3
    SAM_CIGAR = 5
    SAM_SEQ = 9
    # Read sam from stdin
    for line in samfile:
        # Skip header rows
        if line[0][0] == '@':
            continue
        samline = line.split('\t')
        # Check that aln contains soft clipping
        if 'S' in samline[SAM_CIGAR] and 'H' not in samline[SAM_CIGAR]:
            # Get L/R clip lengths
            leftClipLen, rightClipLen = checkClips(samline[SAM_CIGAR])
            alnLen = lenCIGAR(samline[SAM_CIGAR])
            # Check for left overhang
            if leftClipLen:
                if (int(samline[SAM_POS]) <= max_break) and (
                    leftClipLen >= (int(samline[SAM_POS]) + min_clip)
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
                                'L',
                            )
                        )
                    except KeyError:
                        logging.warning(
                            'Reference sequence not found in FAI file: '
                            + str(samline[SAM_RNAME])
                        )
            # Check for right overhang
            if rightClipLen:
                try:
                    ContigLen = contigs[str(samline[SAM_RNAME])]
                except KeyError:
                    logging.warning(
                        'Reference sequence not found in FAI file: '
                        + str(samline[SAM_RNAME])
                    )
                alnEnd = int(samline[SAM_POS]) + alnLen
                # Check if overhang is on contig right end
                if ((ContigLen - alnEnd) <= max_break) and (
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
                            'R',
                        )
                    )


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
    outpaths = []
    # Counter
    readCount = 0
    # For each aligned read
    for read in alignments:
        readCount += 1
        # Log update every 10K reads
        if not readCount % 10000:
            logging.info('Alignments processed: %s' % str(readCount))
        # Check if alignment is at right end overhang
        if read[6] == 'R':
            # i.e [(alnStart,alnEnd,rightClipLen,readSeq,readName,contigName,'R')]
            if prefix:
                base = '_'.join([prefix, str(read[5])])
            else:
                base = str(read[5])
            # Set output file path
            outfileR = os.path.join(outdir, '_'.join([base, 'R']) + '.fasta')
            # Check if need to create new outfile or append to existing file.
            if outfileR not in outpaths:
                logging.info('Creating new outfile: %s' % str(outfileR))
                outpaths.append(outfileR)
                # Output reads aligned to right end of contig
                with open(outfileR, 'w') as fileR:
                    # Convert overhanging segment of read to lowercase
                    masked = read[3][: -read[2]] + read[3][-read[2] :].lower()
                    # Write masked read to fasta
                    writefasta(fileR, str(read[4]), masked)
            else:
                # If file already exists
                # Output reads aligned to right end of contig
                with open(outfileR, 'a') as fileR:
                    # Convert overhanging segment of read to lowercase
                    masked = read[3][: -read[2]] + read[3][-read[2] :].lower()
                    # Write masked read to fasta
                    writefasta(fileR, str(read[4]), masked)
        # Check if alignment is a left end overhang
        elif read[6] == 'L':
            # i.e [(alnStart,alnEnd,leftClipLen,readSeq,readname,contigname,'L')]
            if prefix:
                base = '_'.join([prefix, str(read[5])])
            else:
                base = str(read[5])
            # Set output file path
            outfileL = os.path.join(outdir, '_'.join([base, 'L']) + '.fasta')
            # Check if need to create new outfile or append to existing file.
            if outfileL not in outpaths:
                logging.info('Creating new outfile: %s' % str(outfileL))
                outpaths.append(outfileL)
                # Output reads aligned to left end of contig
                with open(outfileL, 'w') as fileL:
                    # Convert overhanging segment of read to lowercase
                    masked = read[3][: read[2]].lower() + read[3][read[2] :]
                    # Write masked read to fasta
                    writefasta(fileL, str(read[4]), masked)
            else:
                # If file already exists
                # Output reads aligned to left end of contig
                with open(outfileL, 'a') as fileL:
                    # Convert overhanging segment of read to lowercase
                    masked = read[3][: read[2]].lower() + read[3][read[2] :]
                    # Write masked read to fasta
                    writefasta(fileL, str(read[4]), masked)

    logging.info('Total alignments processed: %s' % str(readCount))


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
