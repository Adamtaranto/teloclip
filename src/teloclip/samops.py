"""
SAM file operations for Teloclip.
"""

import logging
import re
import sys
from typing import TYPE_CHECKING, Dict, Iterator, Optional, Tuple

from teloclip.motifs import check_sequence_for_patterns
from teloclip.seqops import isMotifInClip

if TYPE_CHECKING:
    from .extract_io import ExtractionStats


def processSamlines(
    samfile,
    contig_dict,
    motif_list=None,
    match_anywhere=False,
    max_break=0,
    min_clip=1,
    min_repeats=1,
    min_anchor=100,
    return_counts=False,
    exclude_secondary=True,
):
    """
    Process SAM alignment lines and filter based on clipping and motif criteria.

    This function reads SAM format alignment data and filters alignments based on
    soft-clipping patterns at contig ends, motif matching, and quality thresholds.

    Parameters
    ----------
    samfile : file-like object
        Input SAM format file or stream.
    contig_dict : dict
        Dictionary mapping contig names to their lengths.
    motif_list : list of str, optional
        List of motif patterns to search for in sequences. Default is None.
    match_anywhere : bool, optional
        If True, search for motifs anywhere in the sequence. If False, search
        only in clipped regions. Default is False.
    max_break : int, optional
        Maximum allowed gap in alignment. Default is 0.
    min_clip : int, optional
        Minimum soft-clipping length required. Default is 1.
    min_repeats : int, optional
        Minimum number of motif repeats required for a match. Default is 1.
    min_anchor : int, optional
        Minimum anchored alignment length required. Default is 500.
    return_counts : bool, optional
        If True, return processing statistics as a dictionary instead of None.
        Default is False for backward compatibility.
    exclude_secondary : bool, optional
        If True, exclude secondary alignments (FLAG & 256). Default is True.

    Returns
    -------
    None or dict
        If return_counts is False (default), function processes input and writes
        filtered results to stdout, returning None. If return_counts is True,
        returns a dictionary with processing statistics including samlineCount,
        keepCount, motifCount, removeCount, anchorFilteredCount, and bothCount.
        Logging information is provided about processing statistics.
    """
    # Compile motif regex patterns and include min_repeats
    # for faster repeated matching.
    if motif_list:
        if min_repeats > 1:
            logging.info(f'Applying minimum repeats filter: {min_repeats}')
            motif_list = [
                rf'({motif})' + r'{' + f'{min_repeats},' + r'}' for motif in motif_list
            ]
        compiled_motifs = [re.compile(motif) for motif in motif_list]
    else:
        compiled_motifs = []

    if compiled_motifs:
        logging.info(
            f'Compiled motif patterns: {", ".join([str(motif) for motif in compiled_motifs])}'
        )

    # SAM line index keys
    SAM_QNAME = 0
    SAM_FLAG = 1
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

    # Exclusion criteria counters
    excluded_unmapped = 0
    excluded_secondary = 0
    excluded_min_clip = 0
    excluded_max_break = 0
    excluded_min_anchor = 0
    excluded_motifs = 0

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

        # Check for unmapped reads (FLAG & 4)
        flag = int(samline[SAM_FLAG])
        if flag & 4:
            excluded_unmapped += 1
            removeCount += 1
            continue

        # Check for secondary alignments (FLAG & 256)
        if flag & 256 and exclude_secondary:
            excluded_secondary += 1
            removeCount += 1
            continue

        # Check if line contains soft-clip and no hard-clipping.
        if 'S' in samline[SAM_CIGAR] and 'H' not in samline[SAM_CIGAR]:
            # Check if alignment meets minimum anchor requirement
            if not validate_min_anchor(samline[SAM_CIGAR], min_anchor):
                excluded_min_anchor += 1
                anchorFilteredCount += 1
                removeCount += 1
                continue

            # Get length of left and right overhangs
            leftClipLen, rightClipLen = checkClips(samline[SAM_CIGAR])
            alnLen = lenCIGAR(samline[SAM_CIGAR])

            # Track exclusion reasons for reads with clips
            left_excluded_max_break = False
            left_excluded_min_clip = False
            right_excluded_max_break = False
            right_excluded_min_clip = False

            # Check for left overhang
            if leftClipLen:
                pos = int(samline[SAM_POS])
                if pos > max_break:
                    left_excluded_max_break = True
                elif leftClipLen < (pos + min_clip):
                    left_excluded_min_clip = True
                else:
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
                if (ContigLen - alnEnd) > max_break:
                    right_excluded_max_break = True
                elif alnEnd + rightClipLen < ContigLen + 1:
                    right_excluded_min_clip = True
                else:
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

            # Track exclusion reasons for reads that weren't kept
            if not keepLine and (leftClipLen or rightClipLen):
                if left_excluded_max_break or right_excluded_max_break:
                    excluded_max_break += 1
                elif left_excluded_min_clip or right_excluded_min_clip:
                    excluded_min_clip += 1

            # Optional check for Telomeric repeat motifs
            if motif_list and keepLine and match_anywhere:
                if check_sequence_for_patterns(
                    samline[SAM_SEQ], motif_list, min_repeats
                ):
                    sys.stdout.write(line)
                    motifCount += 1
                else:
                    excluded_motifs += 1
                    removeCount += 1
            elif motif_list and keepLine:
                if isMotifInClip(
                    samline,
                    compiled_motifs,
                    leftClip,
                    rightClip,
                    leftClipLen,
                    rightClipLen,
                    1,  # min_repeats already applied in compiled_motifs
                ):
                    sys.stdout.write(line)
                    motifCount += 1
                else:
                    excluded_motifs += 1
                    removeCount += 1
            elif keepLine:
                sys.stdout.write(line)
            else:
                removeCount += 1
    if motif_list:
        logging.info(
            f'Processed {samlineCount} SAM records.\n'
            f'Found {keepCount} alignments soft-clipped at contig ends.\n'
            f'Found {bothCount} alignments spanning entire contigs.\n'
            f'Output {motifCount} alignments containing motif matches.\n'
            f'Exclusion summary:\n'
            f'  - Unmapped reads: {excluded_unmapped}\n'
            f'  - Secondary alignments: {excluded_secondary}\n'
            f'  - Below min_anchor threshold ({min_anchor}bp): {excluded_min_anchor}\n'
            f'  - Beyond max_break threshold ({max_break}bp): {excluded_max_break}\n'
            f'  - Below min_clip threshold: {excluded_min_clip}\n'
            f'  - No telomeric motifs: {excluded_motifs}\n'
            f'Total discarded: {removeCount} alignments after all filtering.'
        )
    else:
        logging.info(
            f'Processed {samlineCount} SAM records.\n'
            f'Found {keepCount} alignments soft-clipped at contig ends.\n'
            f'Found {bothCount} alignments spanning entire contigs.\n'
            f'Exclusion summary:\n'
            f'  - Unmapped reads: {excluded_unmapped}\n'
            f'  - Secondary alignments: {excluded_secondary}\n'
            f'  - Below min_anchor threshold ({min_anchor}bp): {excluded_min_anchor}\n'
            f'  - Beyond max_break threshold ({max_break}bp): {excluded_max_break}\n'
            f'  - Below min_clip threshold: {excluded_min_clip}\n'
            f'Total discarded: {removeCount} alignments after all filtering.'
        )

    # Return counts if requested for testing purposes
    if return_counts:
        return {
            'samlineCount': samlineCount,
            'keepCount': keepCount,
            'motifCount': motifCount,
            'removeCount': removeCount,
            'anchorFilteredCount': anchorFilteredCount,
            'bothCount': bothCount,
            'excluded_unmapped': excluded_unmapped,
            'excluded_secondary': excluded_secondary,
            'excluded_min_anchor': excluded_min_anchor,
            'excluded_max_break': excluded_max_break,
            'excluded_min_clip': excluded_min_clip,
            'excluded_motifs': excluded_motifs,
        }


def splitCIGAR(SAM_CIGAR):
    """
    Split CIGAR string into list of tuples with format (length, operator).

    Parses a CIGAR string and converts it into a list of tuples where each
    tuple contains the length and operation type for each CIGAR element.

    Parameters
    ----------
    SAM_CIGAR : str
        CIGAR string from SAM alignment format.

    Returns
    -------
    list of tuple
        List of tuples with format (length, operator), where length is int
        and operator is str. For example, '174M76S' becomes [(174, 'M'), (76, 'S')].

    Examples
    --------
    >>> splitCIGAR('174M76S')
    [(174, 'M'), (76, 'S')]
    >>> splitCIGAR('96S154M')
    [(96, 'S'), (154, 'M')]
    """
    CIGARlist = []
    for x in re.findall('[0-9]*[A-Z|=]', SAM_CIGAR):
        CIGARlist.append((int(re.findall('[0-9]*', x)[0]), re.findall('[A-Z]|=', x)[0]))
    # 174M76S --> [(174,M),(76,S)]
    # 96S154M --> [(96,S),(154,M)]
    return CIGARlist


def checkClips(SAM_CIGAR):
    """
    Get lengths of soft-clipped blocks from either end of an alignment.

    Analyzes a CIGAR string to determine the lengths of soft-clipped regions
    at the start and end of the alignment.

    Parameters
    ----------
    SAM_CIGAR : str
        CIGAR string from SAM alignment format.

    Returns
    -------
    tuple of (int or None, int or None)
        Tuple containing (left_clip_length, right_clip_length). Returns None
        for positions where no soft-clipping is present.

    Examples
    --------
    >>> checkClips('10S100M20S')
    (10, 20)
    >>> checkClips('100M')
    (None, None)
    >>> checkClips('50S100M')
    (50, None)
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
    Calculate the length of alignment in reference sequence.

    Calculates the total length of the alignment on the reference sequence by
    summing match, deletion, splice, mismatch, and sequence match block values.
    Ignores insertions, padding, hard clips, and soft clips.

    Parameters
    ----------
    SAM_CIGAR : str
        CIGAR string from SAM alignment format.

    Returns
    -------
    int
        Total alignment length on the reference sequence in base pairs.

    Notes
    -----
    Includes operations: M (match), D (deletion), N (splice), X (mismatch), = (sequence match)
    Excludes operations: I (insertion), P (padding), H (hard clip), S (soft clip)
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

    Parameters
    ----------
    cigar_string : str
        CIGAR string from SAM alignment format.

    Returns
    -------
    int
        Number of aligned bases (M + = + X operations only).
    """
    aligned_bases = 0
    cigar_list = splitCIGAR(cigar_string)
    for length, operation in cigar_list:
        if operation in {'M', '=', 'X'}:
            aligned_bases += length
    return aligned_bases


def validate_min_anchor(cigar_string, min_anchor):
    """
    Validate that an alignment has sufficient anchored bases to meet requirement.

    Checks if the alignment has enough actually aligned/matched bases to meet
    the minimum anchor threshold for quality filtering.

    Parameters
    ----------
    cigar_string : str
        CIGAR string from SAM alignment format.
    min_anchor : int
        Minimum number of aligned bases required.

    Returns
    -------
    bool
        True if alignment meets min_anchor requirement, False otherwise.
    """
    aligned_bases = calculate_aligned_bases(cigar_string)
    return aligned_bases >= min_anchor


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


class EnhancedStreamingSamFilter:
    """
    Enhanced streaming SAM filter with additional validation and statistics.

    Improvements over original StreamingSamFilter:
    - Better error handling and validation
    - Statistics tracking
    - Quality and anchor filtering
    - Motif analysis integration

    Parameters
    ----------
    samfile : file-like
        SAM file handle or iterator.
    contigs : Dict[str, int]
        Dictionary of contig names to lengths.
    max_break : int, optional
        Maximum gap from contig end to allow. Default is 50.
    min_clip : int, optional
        Minimum clip length required. Default is 1.
    min_anchor : int, optional
        Minimum anchored alignment length required. Default is 500.
    min_mapq : int, optional
        Minimum mapping quality required. Default is 0.
    motif_patterns : Dict[str, str], optional
        Compiled motif regex patterns. Default is None.
    stats : ExtractionStats, optional
        Statistics tracker. Default is None.
    exclude_secondary : bool, optional
        If True, exclude secondary alignments. Default is True.
    """

    def __init__(
        self,
        samfile,
        contigs: Dict[str, int],
        max_break: int = 50,
        min_clip: int = 1,
        min_anchor: int = 500,
        min_mapq: int = 0,
        motif_patterns: Optional[Dict[str, str]] = None,
        stats: Optional['ExtractionStats'] = None,
        exclude_secondary: bool = True,
    ):
        """
        Initialize enhanced streaming filter.

        Parameters
        ----------
        samfile : file-like
            SAM file handle or iterator
        contigs : Dict[str, int]
            Dictionary of contig names to lengths
        max_break : int
            Maximum gap from contig end to allow
        min_clip : int
            Minimum clip length required
        min_anchor : int
            Minimum anchored alignment length required
        min_mapq : int
            Minimum mapping quality required
        motif_patterns : Dict[str, str], optional
            Compiled motif regex patterns
        stats : ExtractionStats, optional
            Statistics tracker
        exclude_secondary : bool, optional
            If True, exclude secondary alignments. Default is True.
        """
        self.samfile = samfile
        self.contigs = contigs
        self.max_break = max_break
        self.min_clip = min_clip
        self.min_anchor = min_anchor
        self.min_mapq = min_mapq
        self.motif_patterns = motif_patterns or {}
        self.stats = stats or None

        # SAM field indices
        self.SAM_QNAME = 0
        self.SAM_FLAG = 1
        self.SAM_RNAME = 2
        self.SAM_POS = 3
        self.SAM_MAPQ = 4
        self.SAM_CIGAR = 5
        self.SAM_SEQ = 9

    def _split_cigar(self, cigar_string: str) -> list:
        """
        Split CIGAR string into list of (length, operation) tuples.

        Parameters
        ----------
        cigar_string : str
            CIGAR string from SAM alignment.

        Returns
        -------
        list
            List of (length, operation) tuples.
        """
        cigar_list = []
        for match in re.findall(r'[0-9]*[A-Z=]', cigar_string):
            length = int(re.findall(r'[0-9]*', match)[0])
            operation = re.findall(r'[A-Z=]', match)[0]
            cigar_list.append((length, operation))
        return cigar_list

    def _check_clips(self, cigar_string: str) -> Tuple[Optional[int], Optional[int]]:
        """
        Get lengths of soft-clipped blocks from either end of alignment.

        Parameters
        ----------
        cigar_string : str
            CIGAR string from SAM alignment.

        Returns
        -------
        Tuple[Optional[int], Optional[int]]
            Left and right clip lengths (None if no clipping).
        """
        left_clip_len = None
        right_clip_len = None
        cigar_list = self._split_cigar(cigar_string)

        # Check if first segment is soft-clipped
        if cigar_list and cigar_list[0][1] == 'S':
            left_clip_len = cigar_list[0][0]

        # Check if last segment is soft-clipped
        if cigar_list and cigar_list[-1][1] == 'S':
            right_clip_len = cigar_list[-1][0]

        return (left_clip_len, right_clip_len)

    def _calculate_alignment_length(self, cigar_string: str) -> int:
        """
        Calculate alignment length on reference sequence.

        Parameters
        ----------
        cigar_string : str
            CIGAR string from SAM alignment.

        Returns
        -------
        int
            Total alignment length on reference sequence.
        """
        aln_len = 0
        cigar_list = self._split_cigar(cigar_string)
        for length, operation in cigar_list:
            if operation in {'D', 'M', 'N', 'X', '='}:
                aln_len += length
        return aln_len

    def _count_motifs_in_sequence(self, sequence: str) -> Dict[str, int]:
        """
        Count motif occurrences in sequence.

        Parameters
        ----------
        sequence : str
            DNA sequence to search for motifs.

        Returns
        -------
        Dict[str, int]
            Dictionary mapping motif names to occurrence counts.
        """
        motif_counts = {}
        for motif_name, pattern in self.motif_patterns.items():
            matches = re.findall(pattern, sequence)
            motif_counts[motif_name] = len(matches)
        return motif_counts

    def __iter__(self):
        """
        Iterate through filtered SAM alignments.

        Yields
        ------
        dict
            Dictionary containing alignment information with keys:
            - aln_start: int, alignment start position
            - aln_end: int, alignment end position
            - clip_length: int, soft-clip length
            - sequence: str, read sequence
            - read_name: str, read identifier
            - contig_name: str, reference contig name
            - end: str, overhang direction ('L' or 'R')
            - mapq: int, mapping quality score
            - motif_counts: dict, motif occurrence counts (if patterns provided)
            - overhang_seq: str, overhanging sequence portion
        """
        for line in self.samfile:
            # Skip header rows
            if line.startswith('@'):
                continue

            # Count all non-header SAM lines processed
            if self.stats:
                self.stats.record_sam_line()

            try:
                samline = line.strip().split('\t')

                # Basic validation
                if len(samline) < 11:
                    logging.warning(f'Malformed SAM line: {line.strip()}')
                    if self.stats:
                        self.stats.record_filter('malformed')
                    continue

                # Check for unmapped reads (FLAG & 4)
                flag = int(samline[self.SAM_FLAG])
                if flag & 4:
                    if self.stats:
                        self.stats.record_filter('unmapped')
                    continue

                # Check for secondary alignments (FLAG & 256)
                if flag & 256 and self.exclude_secondary:
                    if self.stats:
                        self.stats.record_filter('secondary')
                    continue

                # Check for soft clipping (and no hard clipping)
                cigar = samline[self.SAM_CIGAR]
                if 'S' not in cigar or 'H' in cigar:
                    if self.stats:
                        self.stats.record_filter('soft_clip')
                    continue

                # Quality filtering
                try:
                    mapq = int(samline[self.SAM_MAPQ])
                    if mapq < self.min_mapq:
                        if self.stats:
                            self.stats.record_filter('quality')
                        continue
                except (ValueError, IndexError):
                    logging.warning(f'Invalid MAPQ in line: {line.strip()}')
                    continue

                # Anchor length validation
                if not validate_min_anchor(cigar, self.min_anchor):
                    if self.stats:
                        self.stats.record_filter('anchor')
                    continue

                # Get clip lengths and alignment info
                left_clip_len, right_clip_len = self._check_clips(cigar)
                aln_len = self._calculate_alignment_length(cigar)

                contig_name = samline[self.SAM_RNAME]
                if contig_name not in self.contigs:
                    logging.warning(f'Unknown contig in SAM: {contig_name}')
                    continue

                contig_len = self.contigs[contig_name]
                pos = int(samline[self.SAM_POS])
                sequence = samline[self.SAM_SEQ]
                read_name = samline[self.SAM_QNAME]

                # Track exclusion reasons and process overhangs

                # Check for left overhang
                if left_clip_len:
                    if left_clip_len < self.min_clip:
                        # Track min_clip exclusion but continue to check right overhang
                        if self.stats:
                            self.stats.record_filter('min_clip')
                    elif pos > self.max_break:
                        # Track max_break exclusion but continue to check right overhang
                        if self.stats:
                            self.stats.record_filter('max_break')
                    elif left_clip_len >= (pos + self.min_clip):
                        # Valid left overhang
                        aln_end = pos + aln_len
                        # Extract left overhang sequence (clipped region) and convert to uppercase for motif counting
                        overhang_seq = sequence[:left_clip_len].upper()

                        # Count motifs in the clipped overhang region only
                        motif_counts = None
                        if self.motif_patterns:
                            motif_counts = self._count_motifs_in_sequence(overhang_seq)
                            # Check if no motifs found and track exclusion
                            if not any(count > 0 for count in motif_counts.values()):
                                if self.stats:
                                    self.stats.record_filter('motifs')
                                continue  # Skip this overhang

                        yield {
                            'aln_start': pos,
                            'aln_end': aln_end,
                            'clip_length': left_clip_len,
                            'sequence': sequence,
                            'read_name': read_name,
                            'contig_name': contig_name,
                            'end': 'L',
                            'mapq': mapq,
                            'motif_counts': motif_counts,
                            'overhang_seq': overhang_seq,
                        }
                    else:
                        # Does not meet min_clip requirement relative to position
                        if self.stats:
                            self.stats.record_filter('min_clip')

                # Check for right overhang
                if right_clip_len:
                    aln_end = pos + aln_len
                    if right_clip_len < self.min_clip:
                        # Track min_clip exclusion
                        if self.stats:
                            self.stats.record_filter('min_clip')
                    elif (contig_len - aln_end) > self.max_break:
                        # Track max_break exclusion
                        if self.stats:
                            self.stats.record_filter('max_break')
                    elif aln_end + right_clip_len >= contig_len + 1:
                        # Valid right overhang
                        # Extract right overhang sequence (clipped region) and convert to uppercase for motif counting
                        overhang_seq = sequence[-right_clip_len:].upper()

                        # Count motifs in the clipped overhang region only
                        motif_counts = None
                        if self.motif_patterns:
                            motif_counts = self._count_motifs_in_sequence(overhang_seq)
                            # Check if no motifs found and track exclusion
                            if not any(count > 0 for count in motif_counts.values()):
                                if self.stats:
                                    self.stats.record_filter('motifs')
                                continue  # Skip this overhang

                        yield {
                            'aln_start': pos,
                            'aln_end': aln_end,
                            'clip_length': right_clip_len,
                            'sequence': sequence,
                            'read_name': read_name,
                            'contig_name': contig_name,
                            'end': 'R',
                            'mapq': mapq,
                            'motif_counts': motif_counts,
                            'overhang_seq': overhang_seq,
                        }

            except (IndexError, ValueError) as e:
                logging.warning(f'Error processing SAM line: {e}')
                continue


def enhanced_streaming_split_by_contig(
    alignments: Iterator,
    output_dir: Optional[str] = None,
    prefix: Optional[str] = None,
    output_format: str = 'fasta',
    buffer_size: int = 1000,
    include_stats: bool = False,
    mask_overhangs: bool = True,
    existing_stats: Optional['ExtractionStats'] = None,
    use_sam_attributes: bool = False,
) -> 'ExtractionStats':
    """
    Efficiently write overhang reads using file handles and buffering.

    This is a complete rewrite of StreamingSplitByContig with major improvements:
    - Uses file handle management instead of opening/closing files repeatedly
    - BioPython integration for reliable FASTA/FASTQ writing
    - Buffered writes for better I/O performance
    - Rich sequence headers with optional statistics
    - Motif analysis integration
    - Comprehensive statistics tracking

    Parameters
    ----------
    alignments : Iterator
        Iterator of alignment dictionaries from EnhancedStreamingSamFilter.
    output_dir : str, optional
        Output directory for files.
    prefix : str, optional
        Prefix for output filenames.
    output_format : str
        Output format ('fasta' or 'fastq').
    buffer_size : int
        Number of sequences to buffer before writing.
    include_stats : bool
        Whether to include statistics in sequence headers.
    mask_overhangs : bool
        Whether to convert overhang sequences to lowercase.
    existing_stats : ExtractionStats, optional
        Existing statistics object to update (preserves SAM line counts).
    use_sam_attributes : bool
        Whether to format statistics as SAM attributes for FASTQ output.

    Returns
    -------
    ExtractionStats
        Statistics about the extraction process.
    """
    from .extract_io import ExtractionStats, MultiFileSequenceWriter

    # Use existing stats if provided, otherwise create new one
    stats = existing_stats if existing_stats is not None else ExtractionStats()

    # Initialize multi-file writer
    with MultiFileSequenceWriter(
        base_dir=output_dir,
        prefix=prefix,
        output_format=output_format,
        buffer_size=buffer_size,
        use_sam_attributes=use_sam_attributes,
    ) as writer:
        read_count = 0

        for alignment in alignments:
            read_count += 1

            # Log progress
            if read_count % 10000 == 0:
                logging.info(f'Alignments processed: {read_count}')

            contig_name = alignment['contig_name']
            end = alignment['end']
            sequence = alignment['sequence']
            clip_length = alignment['clip_length']

            # Apply masking if requested
            if mask_overhangs:
                if end == 'L':
                    # Mask left overhang (lowercase)
                    masked_seq = sequence[:clip_length].lower() + sequence[clip_length:]
                else:  # end == 'R'
                    # Mask right overhang (lowercase)
                    masked_seq = (
                        sequence[:-clip_length] + sequence[-clip_length:].lower()
                    )
            else:
                masked_seq = sequence

            # Build sequence description
            description = f'overhang_{end}_{contig_name}'

            # Prepare statistics for header
            seq_stats = None
            if include_stats:
                seq_stats = {
                    'mapq': alignment['mapq'],
                    'clip_length': clip_length,
                    'overhang_length': len(alignment['overhang_seq']),
                }
                if alignment.get('motif_counts'):
                    seq_stats['motif_counts'] = alignment['motif_counts']

            # Write sequence
            writer.write_sequence(
                contig_name=contig_name,
                end=end,
                seq_id=alignment['read_name'],
                sequence=masked_seq,
                description=description,
                stats=seq_stats,
            )

            # Update statistics
            stats.record_alignment(
                contig_name=contig_name,
                is_left=(end == 'L'),
                motif_counts=alignment.get('motif_counts'),
            )

    logging.info(f'Total alignments processed: {stats.total_sam_lines}')
    logging.info(f'Alignments with valid overhangs: {read_count}')
    return stats
