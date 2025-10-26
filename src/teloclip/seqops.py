"""
Sequence operations and utilities for teloclip.

This module provides functions for sequence manipulation, file I/O operations,
and motif analysis in DNA sequences. Includes utilities for FASTA processing,
sequence transformations, and clipping analysis.
"""

from teloclip.motifs import check_sequence_for_patterns
from teloclip.utils import isfile


def makeMask(killIdx, listlen):
    """
    Create a binary mask list with specified indices set to 0.

    Parameters
    ----------
    killIdx : list of int
        List of indices to set to 0 in the mask.
    listlen : int
        Length of the mask list to create.

    Returns
    -------
    list of int
        Binary mask list where specified indices are 0 and others are 1.

    Examples
    --------
    >>> makeMask([0,9], 10)
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 0]
    """
    mask = [1 for i in range(listlen)]
    for x in killIdx:
        mask[x] = 0
    return mask


def filterList(data, exclude):
    """
    Filter a list by excluding elements at specified indices.

    Parameters
    ----------
    data : list
        Input data list to filter.
    exclude : list of int
        List of indices to exclude from the output.

    Returns
    -------
    generator
        Generator yielding elements from data excluding those at specified indices.

    Examples
    --------
    >>> list(filterList([1,2,3,4,5,6,7,8,9,10], [0,9]))
    [2, 3, 4, 5, 6, 7, 8, 9]
    """
    mask = makeMask(exclude, len(data))
    return (d for d, s in zip(data, mask) if s)


def revComp(seq):
    """
    Generate reverse complement of a DNA sequence.

    Parameters
    ----------
    seq : str
        Input DNA sequence string.

    Returns
    -------
    str
        Reverse complement of the input DNA sequence.

    Examples
    --------
    >>> revComp('ATCG')
    'CGAT'
    >>> revComp('AAATTTCCCGGG')
    'CCCGGGAAATTT'
    """

    def revcompl(x):
        """
        Generate reverse complement of DNA sequence.

        Parameters
        ----------
        x : str
            Input DNA sequence string.

        Returns
        -------
        str
            Reverse complement of input sequence.
        """
        return ''.join(
            [{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}[B] for B in x][::-1]
        )

    return revcompl(seq)


def read_fai(fai):
    """
    Import FASTA index file and return dictionary of sequence names and lengths.

    Parameters
    ----------
    fai : str
        Path to FASTA index (.fai) file.

    Returns
    -------
    dict
        Dictionary where keys are sequence names (str) and values are
        sequence lengths (int).
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
    Create unique set of DNA motif strings and their reverse complements.

    Parameters
    ----------
    motifList : list of str
        List of DNA motif strings.

    Returns
    -------
    set of str
        Unique set containing all input motifs and their reverse complements.

    Examples
    --------
    >>> sorted(addRevComplement(['ATCG']))
    ['ATCG', 'CGAT']
    """

    def revcompl(x):
        """
        Generate reverse complement of DNA sequence.

        Parameters
        ----------
        x : str
            Input DNA sequence string.

        Returns
        -------
        str
            Reverse complement of input sequence.
        """
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
    Test for presence of DNA motifs in soft-clipped regions of a read.

    Extracts terminal soft-clipped blocks from read sequence and searches
    for any DNA motif patterns from the provided motif list.

    Parameters
    ----------
    samline : list
        SAM format alignment line split into fields.
    motifList : list of str
        List of DNA motif patterns to search for.
    leftClip : bool
        Whether left clipping is present.
    rightClip : bool
        Whether right clipping is present.
    leftClipLen : int or None
        Length of left soft-clipped region.
    rightClipLen : int or None
        Length of right soft-clipped region.
    minRepeats : int, optional
        Minimum number of motif repeats required for a match. Default is 1.

    Returns
    -------
    bool
        True if any clipped end sequence contains at least one instance
        of any motif, False otherwise.
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
