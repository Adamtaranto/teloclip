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
    revcompl = lambda x: "".join(
        [{"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}[B] for B in x][::-1]
    )
    return revcompl(seq)


# NCU
def writeClip(idx, zpad, gap, seq, maplen):
    # leftpad idx ID
    padIdx = str(idx).zfill(zpad) + ":"
    # If gap between aln end(R) or start(L) and contig end, left pad softclip with '-'
    padseq = "-" * gap + seq
    # Format length of ref covered by alingment
    readlen = "LEN=" + str(maplen).rjust(6)
    print("\t".join([padIdx, readlen, padseq]))


# NCU
def fasta2dict(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    contigDict = dict()
    for header in faiter:
        # Drop the ">"
        # Split on whitespace and take first item as name
        header = header.__next__()[1:].strip()
        name = header.split()[0]
        # Join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        contigDict[name] = (header, seq)
    return contigDict


def writefasta(outfile, name, seq, length=80):
    outfile.write(">" + str(name) + "\n")
    while len(seq) > 0:
        outfile.write(seq[:length] + "\n")
        seq = seq[length:]


# NCU
def manageTemp(record=None, tempPath=None, scrub=False):
    """Create single sequence fasta files or scrub temp files."""
    if scrub and tempPath:
        try:
            os.remove(tempPath)
        except OSError:
            pass
    else:
        with open(tempPath, "w") as f:
            name, seq = record
            writefasta(f, name, seq, length=80)


def read_fai(fai):
    """
    Import fasta index file. Return dict of sequence names and lengths.
    """
    path = isfile(fai)
    # Init empty dict
    ContigDict = dict()
    # Read fai_file to dict
    with open(path, "r") as f:
        for line in f.readlines():
            li = line.strip().split()
            ContigDict[li[0]] = int(li[1])
    return ContigDict


def addRevComplement(motifList):
    """
    Take list of DNA motif strings and return unique set of strings and their reverse complements.
    """
    revcompl = lambda x: "".join(
        [{"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}[B] for B in x][::-1]
    )
    setList = list()
    for motif in motifList:
        setList.append(motif)
        setList.append(revcompl(motif))
    return set(setList)


# Depreciated
def crunchHomopolymers(motifList):
    """
    Take as input a list of target motifs, collapse poly-nucleotide tracks, return list of collapsed motifs.
    """
    # List to catch all collapsed motifs.
    crunchList = list()
    # For each motif
    for motif in motifList:
        # Create empty list to catch not repeated bases.
        noReps = list()
        # Walk through original motif base-by-base.
        for base in motif:
            # If list of kept bases in empty, add first base.
            if not noReps:
                noReps.append(base)
            # If list exists and base is not same as last, store new base.
            elif base != noReps[-1]:
                noReps.append(base)
        # Convert list to string and store new motif
        crunchList.append("".join(noReps))
    # Convert to set to remove duplicates and return
    return list(set(crunchList))


# TODO: support min pattern matches
def isMotifInClip(samline, motifList, leftClip, rightClip, leftClipLen, rightClipLen):
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
            samline[SAM_SEQ][0:leftClipLen], motifList
        )
    if rightClip:
        rightcheck = check_sequence_for_patterns(
            samline[SAM_SEQ][-rightClipLen:], motifList
        )

    # True if either clipped end sequence contains at least one instance of any motif.
    return any([leftcheck, rightcheck])
