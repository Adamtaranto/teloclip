"""
Sequence operations and utilities for teloclip.

This module provides functions for sequence manipulation, file I/O operations,
and motif analysis in DNA sequences. Includes utilities for FASTA processing,
sequence transformations, and clipping analysis.
"""

from teloclip.core.motifs import check_sequence_for_patterns
from teloclip.utils import isfile


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
