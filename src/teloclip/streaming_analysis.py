"""
Memory-efficient analysis functions for processing large genomes.

This module provides streaming analysis functions that process contigs
individually to minimize memory usage.
"""

from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional

import pysam

from .analysis import ContigStats, OverhangInfo
from .samops import checkClips, splitCIGAR


def collect_contig_overhangs_streaming(
    bam_file: pysam.AlignmentFile,
    contig_name: str,
    contig_length: int,
    max_break: int = 10,
    min_clip: int = 1,
    min_anchor: int = 500,
) -> ContigStats:
    """
    Collect overhang statistics for a single contig using streaming access.

    Parameters
    ----------
    bam_file : pysam.AlignmentFile
        Opened BAM file with index.
    contig_name : str
        Name of the contig to analyze.
    contig_length : int
        Length of the contig.
    max_break : int, optional
        Maximum gap allowed between alignment and contig end (default: 10).
    min_clip : int, optional
        Minimum clip length required (default: 1).
    min_anchor : int, optional
        Minimum anchor length required for alignment (default: 500).

    Returns
    -------
    ContigStats
        Statistics for this contig's overhangs.
    """
    contig_stats = ContigStats(contig_name, contig_length)

    try:
        # Fetch alignments only for this contig
        alignments = bam_file.fetch(contig_name)
    except ValueError:
        # Contig not found in BAM file
        return contig_stats

    for alignment in alignments:
        # Skip unmapped reads
        if alignment.is_unmapped:
            continue

        # Skip secondary/supplementary alignments
        if alignment.is_secondary or alignment.is_supplementary:
            continue

        # Check for soft clipping
        cigar_string = alignment.cigarstring
        if not cigar_string or 'S' not in cigar_string or 'H' in cigar_string:
            continue

        # Get clip lengths
        left_clip, right_clip = checkClips(cigar_string)
        left_clip = left_clip or 0
        right_clip = right_clip or 0

        if left_clip == 0 and right_clip == 0:
            continue

        # Calculate alignment positions
        alignment_pos = alignment.reference_start + 1  # Convert to 1-based
        cigar_ops = splitCIGAR(cigar_string)
        alignment_end = (
            alignment_pos + sum(length for length, op in cigar_ops if op in 'MDN=X') - 1
        )

        # Calculate anchor length (aligned portion)
        read_length = alignment.query_length or len(alignment.query_sequence)
        anchor_length = read_length - left_clip - right_clip

        # Check minimum anchor requirement
        if anchor_length < min_anchor:
            continue

        # Check for left clip at contig start
        if left_clip >= min_clip and alignment_pos <= max_break:
            overhang_seq = (
                alignment.query_sequence[:left_clip] if alignment.query_sequence else ''
            )

            overhang = OverhangInfo(
                sequence=overhang_seq,
                length=left_clip,
                alignment_pos=alignment_pos,
                alignment_end=alignment_end,
                read_name=alignment.query_name,
                is_left=True,
                clip_length=left_clip,
                anchor_length=anchor_length,
            )
            contig_stats.left_overhangs.append(overhang)

        # Check for right clip at contig end
        if right_clip >= min_clip and alignment_end >= contig_length - max_break:
            overhang_seq = (
                alignment.query_sequence[-right_clip:]
                if alignment.query_sequence
                else ''
            )

            overhang = OverhangInfo(
                sequence=overhang_seq,
                length=right_clip,
                alignment_pos=alignment_pos,
                alignment_end=alignment_end,
                read_name=alignment.query_name,
                is_left=False,
                clip_length=right_clip,
                anchor_length=anchor_length,
            )
            contig_stats.right_overhangs.append(overhang)

    return contig_stats


def stream_contigs_for_extension(
    bam_file: pysam.AlignmentFile,
    contig_dict: Dict[str, int],
    min_overhangs: int = 1,
    max_break: int = 10,
    min_clip: int = 1,
    min_anchor: int = 500,
    exclude_outliers: bool = False,
    outlier_threshold: float = 2.0,
) -> Iterator[tuple]:
    """
    Stream contigs that meet criteria for extension.

    Parameters
    ----------
    bam_file : pysam.AlignmentFile
        Opened BAM file with index.
    contig_dict : Dict[str, int]
        Dictionary mapping contig names to lengths.
    min_overhangs : int, optional
        Minimum number of overhangs required (default: 1).
    max_break : int, optional
        Maximum gap allowed between alignment and contig end (default: 10).
    min_clip : int, optional
        Minimum clip length required (default: 1).
    min_anchor : int, optional
        Minimum anchor length required (default: 500).
    exclude_outliers : bool, optional
        Whether to exclude outlier contigs (default: False).
    outlier_threshold : float, optional
        Z-score threshold for outlier detection (default: 2.0).

    Yields
    ------
    tuple
        (contig_name, contig_stats) for contigs that meet extension criteria.
    """
    # First pass: collect all stats if outlier detection is needed
    all_stats = {}
    if exclude_outliers:
        for contig_name, contig_length in contig_dict.items():
            stats = collect_contig_overhangs_streaming(
                bam_file, contig_name, contig_length, max_break, min_clip, min_anchor
            )
            all_stats[contig_name] = stats

        # Import here to avoid circular imports
        from .analysis import identify_outlier_contigs

        outliers = identify_outlier_contigs(all_stats, outlier_threshold)
        outlier_set = set(outliers['left_outliers'] + outliers['right_outliers'])

        # Yield non-outlier contigs with sufficient overhangs
        for contig_name, contig_stats in all_stats.items():
            if contig_name not in outlier_set:
                if (
                    len(contig_stats.left_overhangs) >= min_overhangs
                    or len(contig_stats.right_overhangs) >= min_overhangs
                ):
                    yield contig_name, contig_stats
    else:
        # Stream contigs individually without outlier detection
        for contig_name, contig_length in contig_dict.items():
            stats = collect_contig_overhangs_streaming(
                bam_file, contig_name, contig_length, max_break, min_clip, min_anchor
            )

            # Check if this contig has sufficient overhangs for extension
            if (
                len(stats.left_overhangs) >= min_overhangs
                or len(stats.right_overhangs) >= min_overhangs
            ):
                yield contig_name, stats


@dataclass
class ExtensionResult:
    """Result of a contig extension operation."""

    contig_name: str
    original_length: int
    extended_sequence: str
    extension_info: dict
    warnings: List[str] = field(default_factory=list)
    motif_counts: Dict[str, int] = field(default_factory=dict)


def process_single_contig_extension(
    contig_name: str,
    contig_stats: ContigStats,
    original_sequence: str,
    min_extension: int = 1,
    max_homopolymer: int = 50,
    motif_patterns: Optional[Dict[str, str]] = None,
    dry_run: bool = False,
) -> Optional[ExtensionResult]:
    """
    Process extension for a single contig.

    Parameters
    ----------
    contig_name : str
        Name of the contig.
    contig_stats : ContigStats
        Overhang statistics for the contig.
    original_sequence : str
        Original contig sequence.
    min_extension : int, optional
        Minimum extension length required (default: 1).
    max_homopolymer : int, optional
        Maximum homopolymer run allowed (default: 50).
    motif_patterns : Optional[Dict[str, str]], optional
        Motif patterns to search for.
    dry_run : bool, optional
        Whether this is a dry run (default: False).

    Returns
    -------
    Optional[ExtensionResult]
        Extension result if successful, None otherwise.
    """
    # Import here to avoid circular imports
    import re

    from .analysis import detect_homopolymer_runs, select_best_overhang
    from .extension import apply_contig_extension

    warnings = []

    # Check each end for extension opportunities
    for is_left in [True, False]:
        overhangs = (
            contig_stats.left_overhangs if is_left else contig_stats.right_overhangs
        )
        end_name = 'left' if is_left else 'right'

        if not overhangs:
            continue

        # Select best overhang
        best_overhang = select_best_overhang(overhangs, min_extension, max_homopolymer)

        if best_overhang is None:
            continue

        # Check for homopolymer runs
        homo_runs = detect_homopolymer_runs(best_overhang.sequence, max_homopolymer)
        if homo_runs:
            warning_msg = (
                f'Homopolymer run detected in {contig_name} {end_name} extension: '
                f'{homo_runs[0][0]}x{homo_runs[0][3]} at position {homo_runs[0][1]}'
            )
            warnings.append(warning_msg)

        # Apply extension
        if not dry_run:
            try:
                extended_seq, ext_info = apply_contig_extension(
                    original_sequence, best_overhang, contig_stats.contig_length
                )
            except ValueError as e:
                warnings.append(f'Extension failed for {contig_name}: {e}')
                continue
        else:
            # Simulate extension for dry run
            extended_seq = original_sequence  # Don't actually extend
            ext_info = {
                'overhang_length': best_overhang.length,
                'read_name': best_overhang.read_name,
                'is_left': is_left,
                'original_length': contig_stats.contig_length,
                'final_length': contig_stats.contig_length + best_overhang.length,
                'trim_length': 0,
            }

        # Count motifs if patterns provided
        motif_counts = {}
        if motif_patterns:
            target_seq = (
                extended_seq
                if not dry_run
                else original_sequence + best_overhang.sequence
            )
            for pattern_name, pattern_str in motif_patterns.items():
                matches = re.findall(pattern_str, target_seq)
                motif_counts[pattern_name] = len(matches)

        return ExtensionResult(
            contig_name=contig_name,
            original_length=contig_stats.contig_length,
            extended_sequence=extended_seq,
            extension_info=ext_info,
            warnings=warnings,
            motif_counts=motif_counts,
        )

    return None  # No suitable extension found
