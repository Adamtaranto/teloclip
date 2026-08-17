"""SAM file operations for Teloclip."""

from collections import defaultdict
import logging
import re
import sys
from typing import TYPE_CHECKING, Dict, Iterator, Optional

from teloclip.core.motifs import check_sequence_for_patterns
from teloclip.core.overhang import (
    REASON_MAX_BREAK,
    REASON_MIN_CLIP,
    REASON_NO_CLIP,
    classify,
    ends_from_sam_fields,
)
from teloclip.core.seqops import isMotifInClip
from teloclip.io.formats import InputFormatError, iter_sam_lines
from teloclip.report.text import histogram

if TYPE_CHECKING:
    from .extract import ExtractionStats


# Malformed records are reported individually up to this many, then only
# counted. A systematically broken input would otherwise bury every other
# message in the log under one warning per line.
MAX_MALFORMED_WARNINGS = 5


def _warn_malformed(count: int, line: str, reason: str) -> None:
    """
    Log a skipped SAM record, without letting a broken file flood the log.

    Parameters
    ----------
    count : int
        How many malformed records have been seen so far, including this one.
    line : str
        The offending line, truncated in the message.
    reason : str
        Why the line could not be parsed.
    """
    if count > MAX_MALFORMED_WARNINGS:
        if count == MAX_MALFORMED_WARNINGS + 1:
            logging.warning(
                'Further malformed SAM records will be counted but not listed. '
                'The total is reported in the exclusion summary.'
            )
        return

    excerpt = line.rstrip('\n')
    if len(excerpt) > 120:
        excerpt = excerpt[:117] + '...'
    logging.warning(f'Skipping malformed SAM record ({reason}): {excerpt}')


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
    SAM_CIGAR = 5
    SAM_SEQ = 9

    # The SAM specification defines eleven mandatory fields. Anything shorter
    # cannot be indexed safely, whatever else it might be.
    SAM_MIN_FIELDS = 11

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
    excluded_no_clip = 0
    excluded_malformed = 0

    # Per-contig-end tallies of kept overhangs, for the closing summary.
    left_by_contig: Dict[str, int] = defaultdict(int)
    right_by_contig: Dict[str, int] = defaultdict(int)

    # Read SAM from stdin. iter_sam_lines converts a mid-stream decoding
    # failure into a message naming BAM as the likely cause.
    for line in iter_sam_lines(samfile, command='filter'):
        keepLine = False
        leftClip = False
        rightClip = False

        # A blank line is not a record and not an error: a trailing newline at
        # end of file is normal, and neither counting nor warning about it
        # would tell the user anything.
        if not line.strip():
            continue

        # Write headers to stdout
        if line.startswith('@'):
            sys.stdout.write(line)
            continue
        samlineCount += 1
        samline = line.rstrip('\n').split('\t')

        # A truncated or non-SAM line is skipped rather than fatal. A run
        # killed by one bad record part way through a large alignment file
        # loses all the work before it, and a truncated final record is the
        # ordinary result of an interrupted upstream process. The count is
        # reported at the end so a systematically broken input is still
        # obvious rather than silently producing an empty result.
        if len(samline) < SAM_MIN_FIELDS:
            excluded_malformed += 1
            removeCount += 1
            _warn_malformed(excluded_malformed, line, 'too few tab-separated fields')
            continue

        try:
            flag = int(samline[SAM_FLAG])
        except ValueError:
            excluded_malformed += 1
            removeCount += 1
            _warn_malformed(
                excluded_malformed,
                line,
                f'FLAG is not an integer: {samline[SAM_FLAG]!r}',
            )
            continue

        # Check for unmapped reads (FLAG & 4)
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
            try:
                ContigLen = contig_dict[str(samline[SAM_RNAME])]
            except KeyError:
                sys.exit(
                    'Reference sequence not found in FAI file: '
                    + str(samline[SAM_RNAME])
                )

            # POS and CIGAR are only parsed here, so a record with a plausible
            # field count but unparseable coordinates surfaces at this point
            # rather than at the top of the loop.
            try:
                ends = ends_from_sam_fields(samline, ContigLen)
            except (ValueError, IndexError) as error:
                excluded_malformed += 1
                removeCount += 1
                _warn_malformed(excluded_malformed, line, str(error))
                continue

            # Check if alignment meets minimum anchor requirement
            if ends.anchor < min_anchor:
                excluded_min_anchor += 1
                anchorFilteredCount += 1
                removeCount += 1
                continue

            # Both ends are judged by the shared predicate, so filter, extract
            # and extend agree on what counts as a terminal overhang.
            left_call, right_call = classify(
                ends, max_break=max_break, min_clip=min_clip
            )
            leftClipLen = ends.left_clip
            rightClipLen = ends.right_clip

            # Check for left overhang
            if left_call.accepted:
                keepLine = True
                leftClip = True
                keepCount += 1
                left_by_contig[samline[SAM_RNAME]] += 1

            # Check for right overhang
            if right_call.accepted:
                rightClip = True
                right_by_contig[samline[SAM_RNAME]] += 1
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

            # Track exclusion reasons for reads that weren't kept. Every
            # discarded read lands in exactly one bucket, so the buckets sum to
            # the discard total.
            if not keepLine:
                reasons = {left_call.reason, right_call.reason}
                if REASON_MAX_BREAK in reasons:
                    excluded_max_break += 1
                elif REASON_MIN_CLIP in reasons:
                    excluded_min_clip += 1
                else:
                    # Both ends reported no clip: the CIGAR contains an 'S' but
                    # not at either end, so there is no terminal overhang.
                    excluded_no_clip += 1

            # Optional check for Telomeric repeat motifs
            if motif_list and keepLine and match_anywhere:
                # min_repeats is already baked into motif_list as a regex
                # quantifier above, so passing it again here would demand the
                # repeat count twice over.
                if check_sequence_for_patterns(samline[SAM_SEQ], motif_list, 1):
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
        else:
            # No soft clip, or a hard clip whose bases are absent from the
            # record. Previously this fell through with no bucket at all, so
            # the exclusion reasons never summed to the discard total.
            excluded_no_clip += 1
            removeCount += 1

    # One summary, whether or not motifs were requested. The buckets are
    # mutually exclusive and sum to the discard total, so the figures can be
    # checked against each other.
    summary = [
        f'Processed {samlineCount} SAM records.',
        f'Found {keepCount} alignments soft-clipped at contig ends.',
        f'Found {bothCount} alignments spanning entire contigs.',
    ]
    if motif_list:
        summary.append(f'Output {motifCount} alignments containing motif matches.')

    summary.append('Exclusion summary:')
    if excluded_malformed:
        summary.append(f'  - Malformed records (skipped): {excluded_malformed}')
    summary.append(f'  - Unmapped reads: {excluded_unmapped}')
    summary.append(f'  - Secondary alignments: {excluded_secondary}')
    summary.append(f'  - No usable soft clip: {excluded_no_clip}')
    summary.append(
        f'  - Below min_anchor threshold ({min_anchor}bp): {excluded_min_anchor}'
    )
    summary.append(
        f'  - Beyond max_break threshold ({max_break}bp): {excluded_max_break}'
    )
    summary.append(
        f'  - Clip does not reach {min_clip}bp past contig end: {excluded_min_clip}'
    )
    if motif_list:
        summary.append(f'  - No telomeric motifs: {excluded_motifs}')
    summary.append(f'Total discarded: {removeCount} alignments after all filtering.')

    logging.info('\n'.join(summary))

    # Every record unparseable is not a stray bad line, it is the wrong input.
    # Reporting "0 alignments kept" would look like a successful run that found
    # nothing, which is the same output as a clean file with no overhangs.
    if samlineCount and excluded_malformed == samlineCount:
        raise InputFormatError(
            f'None of the {samlineCount} records read could be parsed as SAM. '
            f'Check that the input is uncompressed SAM with the eleven '
            f'mandatory fields, for example the output of `samtools view -h`.'
        )

    # Report the per-contig-end distribution of what survived. A long right
    # tail here is the signature of a collapsed repeat or an organellar contig
    # attracting reads from across the assembly.
    if left_by_contig or right_by_contig:
        per_end_counts = [
            count
            for contig in sorted(set(left_by_contig) | set(right_by_contig))
            for count in (left_by_contig[contig], right_by_contig[contig])
        ]
        for contig in sorted(set(left_by_contig) | set(right_by_contig)):
            logging.info(
                f'{contig}: {left_by_contig[contig]} left, '
                f'{right_by_contig[contig]} right overhang reads'
            )
        chart = histogram(per_end_counts, label='overhang reads')
        if chart:
            logging.info(
                f'Overhang reads per contig end ({len(per_end_counts)} ends):\n{chart}'
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
            'excluded_no_clip': excluded_no_clip,
            'left_by_contig': dict(left_by_contig),
            'right_by_contig': dict(right_by_contig),
        }


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
        """Initialize enhanced streaming filter.

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
        self.exclude_secondary = exclude_secondary

        # SAM field indices
        self.SAM_QNAME = 0
        self.SAM_FLAG = 1
        self.SAM_RNAME = 2
        self.SAM_POS = 3
        self.SAM_MAPQ = 4
        self.SAM_CIGAR = 5
        self.SAM_SEQ = 9

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
        for line in iter_sam_lines(self.samfile, command='extract'):
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

                contig_name = samline[self.SAM_RNAME]
                if contig_name not in self.contigs:
                    logging.warning(f'Unknown contig in SAM: {contig_name}')
                    continue

                ends = ends_from_sam_fields(samline, self.contigs[contig_name])

                # Anchor length validation
                if ends.anchor < self.min_anchor:
                    if self.stats:
                        self.stats.record_filter('anchor')
                    continue

                # Both ends are judged by the shared predicate, so filter,
                # extract and extend agree on what counts as a terminal
                # overhang.
                calls = classify(ends, max_break=self.max_break, min_clip=self.min_clip)

            except (IndexError, ValueError) as e:
                logging.warning(f'Error processing SAM line: {e}')
                continue

            sequence = samline[self.SAM_SEQ]
            read_name = samline[self.SAM_QNAME]

            for call in calls:
                if not call.accepted:
                    # A missing clip is not an exclusion worth reporting: most
                    # reads are clipped at one end only.
                    if call.reason != REASON_NO_CLIP and self.stats:
                        self.stats.record_filter(call.reason)
                    continue

                # Uppercased for motif counting.
                overhang_seq = call.overhang_sequence(ends).upper()

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
                    'aln_start': ends.aln_start,
                    'aln_end': ends.aln_end,
                    'clip_length': call.clip_len,
                    'sequence': sequence,
                    'read_name': read_name,
                    'contig_name': contig_name,
                    'end': 'L' if call.is_left else 'R',
                    'mapq': mapq,
                    'motif_counts': motif_counts,
                    'overhang_seq': overhang_seq,
                }


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
    from .extract import ExtractionStats, MultiFileSequenceWriter

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
