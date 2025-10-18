"""
Enhanced streaming operations for efficient overhang extraction.

This module provides refactored versions of the streaming SAM processing
functions with better performance, error handling, and feature integration.
"""

import logging
import re
from typing import Dict, Iterator, Optional, Tuple

from .extract_io import ExtractionStats, MultiFileSequenceWriter
from .samops import validate_min_anchor


class EnhancedStreamingSamFilter:
    """
    Enhanced streaming SAM filter with additional validation and statistics.

    Improvements over original StreamingSamFilter:
    - Better error handling and validation
    - Statistics tracking
    - Quality and anchor filtering
    - Motif analysis integration
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
        stats: Optional[ExtractionStats] = None,
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
        """
        self.samfile = samfile
        self.contigs = contigs
        self.max_break = max_break
        self.min_clip = min_clip
        self.min_anchor = min_anchor
        self.min_mapq = min_mapq
        self.motif_patterns = motif_patterns or {}
        self.stats = stats or ExtractionStats()

        # SAM field indices
        self.SAM_QNAME = 0
        self.SAM_FLAG = 1
        self.SAM_RNAME = 2
        self.SAM_POS = 3
        self.SAM_MAPQ = 4
        self.SAM_CIGAR = 5
        self.SAM_SEQ = 9

    def _split_cigar(self, cigar_string: str) -> list:
        """Split CIGAR string into list of (length, operation) tuples."""
        cigar_list = []
        for match in re.findall(r'[0-9]*[A-Z=]', cigar_string):
            length = int(re.findall(r'[0-9]*', match)[0])
            operation = re.findall(r'[A-Z=]', match)[0]
            cigar_list.append((length, operation))
        return cigar_list

    def _check_clips(self, cigar_string: str) -> Tuple[Optional[int], Optional[int]]:
        """Get lengths of soft-clipped blocks from either end of alignment."""
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
        """Calculate alignment length on reference sequence."""
        aln_len = 0
        cigar_list = self._split_cigar(cigar_string)
        for length, operation in cigar_list:
            if operation in {'D', 'M', 'N', 'X', '='}:
                aln_len += length
        return aln_len

    def _count_motifs_in_sequence(self, sequence: str) -> Dict[str, int]:
        """Count motif occurrences in sequence."""
        motif_counts = {}
        for motif_name, pattern in self.motif_patterns.items():
            matches = re.findall(pattern, sequence)
            motif_counts[motif_name] = len(matches)
        return motif_counts

    def __iter__(self):
        """Iterate through filtered SAM alignments."""
        for line in self.samfile:
            # Skip header rows
            if line.startswith('@'):
                continue

            try:
                samline = line.strip().split('\\t')

                # Basic validation
                if len(samline) < 11:
                    logging.warning(f'Malformed SAM line: {line.strip()}')
                    continue

                # Check for soft clipping (and no hard clipping)
                cigar = samline[self.SAM_CIGAR]
                if 'S' not in cigar or 'H' in cigar:
                    continue

                # Quality filtering
                try:
                    mapq = int(samline[self.SAM_MAPQ])
                    if mapq < self.min_mapq:
                        self.stats.record_filter('quality')
                        continue
                except (ValueError, IndexError):
                    logging.warning(f'Invalid MAPQ in line: {line.strip()}')
                    continue

                # Anchor length validation
                if not validate_min_anchor(cigar, self.min_anchor):
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

                # Count motifs if patterns provided
                motif_counts = None
                if self.motif_patterns:
                    motif_counts = self._count_motifs_in_sequence(sequence)

                # Check for left overhang
                if left_clip_len and left_clip_len >= self.min_clip:
                    if pos <= self.max_break and left_clip_len >= (pos + self.min_clip):
                        aln_end = pos + aln_len
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
                            'overhang_seq': sequence[:left_clip_len],
                        }

                # Check for right overhang
                if right_clip_len and right_clip_len >= self.min_clip:
                    aln_end = pos + aln_len
                    if (
                        contig_len - aln_end
                    ) <= self.max_break and aln_end + right_clip_len >= contig_len + 1:
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
                            'overhang_seq': sequence[-right_clip_len:],
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
) -> ExtractionStats:
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
        Iterator of alignment dictionaries from EnhancedStreamingSamFilter
    output_dir : str, optional
        Output directory for files
    prefix : str, optional
        Prefix for output filenames
    output_format : str
        Output format ('fasta' or 'fastq')
    buffer_size : int
        Number of sequences to buffer before writing
    include_stats : bool
        Whether to include statistics in sequence headers
    mask_overhangs : bool
        Whether to convert overhang sequences to lowercase

    Returns
    -------
    ExtractionStats
        Statistics about the extraction process
    """
    stats = ExtractionStats()

    # Initialize multi-file writer
    with MultiFileSequenceWriter(
        base_dir=output_dir,
        prefix=prefix,
        output_format=output_format,
        buffer_size=buffer_size,
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

    logging.info(f'Total alignments processed: {read_count}')
    return stats
