from itertools import groupby

from teloclip.motifs import check_sequence_for_patterns
from teloclip.utils import isfile


def makeMask(killIdx, listlen):
    # makeMask([0,9], 10) =  [0,1,1,1,1,1,1,1,1,0]
    mask = [1 for i in range(listlen)]
    for x in killIdx:
        mask[x] = 0
    return mask


# NCU
def filterList(data, exclude):
    # filterList([1,2,3,4,5,6,7,8,9,10],[0,9]) = [2,3,4,5,6,7,8,9]
    mask = makeMask(exclude, len(data))
    return (d for d, s in zip(data, mask) if s)


# NCU
def revComp(seq):
    """Rev comp DNA string."""

    def revcompl(x):
        return ''.join(
            [{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}[B] for B in x][::-1]
        )

    return revcompl(seq)


# NCU
def writeClip(idx, zpad, gap, seq, maplen):
    # leftpad idx ID
    padIdx = str(idx).zfill(zpad) + ':'
    # If gap between aln end(R) or start(L) and contig end, left pad softclip with '-'
    padseq = '-' * gap + seq
    # Format length of ref covered by alingment
    readlen = 'LEN=' + str(maplen).rjust(6)
    print('\t'.join([padIdx, readlen, padseq]))


# NCU
def fasta2dict(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == '>'))
    contigDict = {}
    for header in faiter:
        # Drop the ">"
        # Split on whitespace and take first item as name
        header = header.__next__()[1:].strip()
        name = header.split()[0]
        # Join all sequence lines to one.
        seq = ''.join(s.strip() for s in faiter.__next__())
        contigDict[name] = (header, seq)
    return contigDict


def writefasta(outfile, name, seq, length=80):
    outfile.write('>' + str(name) + '\n')
    while len(seq) > 0:
        outfile.write(seq[:length] + '\n')
        seq = seq[length:]


def read_fai(fai):
    """
    Import fasta index file. Return dict of sequence names and lengths.
    """
    path = isfile(fai)
    # Init empty dict
    ContigDict = {}
    # Read fai_file to dict
    with open(path, 'r') as f:
        for line in f.readlines():
            li = line.strip().split()
            ContigDict[li[0]] = int(li[1])
    return ContigDict


def addRevComplement(motifList):
    """
    Take list of DNA motif strings and return unique set of strings and their reverse complements.
    """

    def revcompl(x):
        return ''.join(
            [{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}[B] for B in x][::-1]
        )

    setList = []
    for motif in motifList:
        setList.append(motif)
        setList.append(revcompl(motif))
    return set(setList)


def isMotifInClip(
    samline, motifList, leftClip, rightClip, leftClipLen, rightClipLen, minRepeats=1
):
    """
    Extract terminal soft-clipped blocks from read sequence and test for presence of any DNA motif in motifList.
    """
    # Sam seq field
    SAM_SEQ = 9

    # Initialize leftcheck and rightcheck
    leftcheck = False
    rightcheck = False

    # Search motif/s as regex in the clipped segment
    if leftClip:
        leftcheck = check_sequence_for_patterns(
            samline[SAM_SEQ][0:leftClipLen], motifList, minRepeats
        )
    if rightClip:
        rightcheck = check_sequence_for_patterns(
            samline[SAM_SEQ][-rightClipLen:], motifList, minRepeats
        )

    # True if either clipped end sequence contains at least one instance of any motif.
    return any([leftcheck, rightcheck])
