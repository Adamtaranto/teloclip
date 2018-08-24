#!/usr/bin/env python
#python 3
#Contact, Adam Taranto, aptaranto@ucdavis.edu

##########################################################################
# Find soft-clipped alignments containing unassembled telomeric repeats. #
##########################################################################

import os
import sys
import re

__version__ = "0.0.2"

class Error (Exception): pass

def isfile(path):
    """
    Test for existence of input file.
    """
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        print("Input file not found: %s" % path)
        sys.exit(1)
    else:
        return path

def read_fai(fai):
    """
    Import fasta index file. Return dict of sequence names and lengths.
    """
    path = isfile(fai)
    # Init empty dict
    ContigDict = dict()
    # Read fai_file to dict
    with open(path, "rU") as f:
        for line in f.readlines():
            li = line.strip().split()
            ContigDict[li[0]]=int(li[1])
    return ContigDict

def addRevComplement(motifList):
    """
    Take list of DNA motif strings and return unique set of strings and their reverse complements.
    """
    revcompl = lambda x: ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}[B] for B in x][::-1])
    setList = list()
    for motif in motifList:
        setList.append(motif)
        setList.append(revcompl(motif))
    return set(setList)

def splitCIGAR(SAM_CIGAR):
    """
    Split CIGAR string into list of tuples with format (len,operator)
    """
    CIGARlist = list()
    for x in re.findall('[0-9]*[A-Z|=]',SAM_CIGAR):
        CIGARlist.append((int(re.findall('[0-9]*', x)[0]),re.findall('[A-Z]|=', x)[0]))
    #174M76S --> [(174,M),(76,S)]
    #96S154M --> [(96,S),(154,M)]
    return CIGARlist

def lenCIGAR(SAM_CIGAR):
    """
    Calculate length of alignment in reference sequence as sum of 
    match, read-deletion, splice, mismatch, and read-match block values.
    Ignore read-insertions, padding, hard and soft clip blocks.
    """
    alnLen = 0
    CIGARlist = splitCIGAR(SAM_CIGAR)
    for x in CIGARlist: # i.e. = [(174,M),(76,S)]
        if x[1] in set(['D','M','N','X','=']):
            alnLen += x[0]
    #Ignore operators in set('P','H','S','I')
    return alnLen

def checkClips(SAM_CIGAR):
    """
    Get lengths of soft-clipped blocks from either end of an alignment given a CIGAR string.
    """
    leftClipLen = None
    rightClipLen = None
    CIGARlist = splitCIGAR(SAM_CIGAR)
    # Check if first segment is soft-clipped
    if CIGARlist[0][1] == "S" :
        leftClipLen = int(CIGARlist[0][0])
    # Check if last segment is soft-clipped
    if CIGARlist[-1][1] == "S" :
        rightClipLen = int(CIGARlist[-1][0])
    return (leftClipLen,rightClipLen)

def isClipMotif(samline,motifList,leftClip,rightClip,leftClipLen,rightClipLen):
    """
    Extract terminal soft-clipped blocks from read sequence and test for presence of any DNA motif in motifList.
    """
    clipSeq = list()
    SAM_SEQ = 9
    if leftClip:
        clipSeq.append(samline[SAM_SEQ][0:leftClipLen])
    if rightClip:
        clipSeq.append(samline[SAM_SEQ][-rightClipLen:])
    # True if either clipped end sequence contains at least one instance of any motif
    return any(s in x for s in motifList for x in clipSeq)

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
