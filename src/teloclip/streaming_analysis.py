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
                contig_name=contig_name,
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
                contig_name=contig_name,
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
    Process extension for a single contig, handling both ends if valid.

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

    from .analysis import select_best_overhang
    from .extension import apply_contig_extension

    warnings = []

    # Find best overhangs for both ends
    best_left_overhang = None
    best_right_overhang = None

    # Check left end
    if contig_stats.left_overhangs:
        best_left_overhang = select_best_overhang(
            contig_stats.left_overhangs, min_extension, max_homopolymer
        )

    # Check right end
    if contig_stats.right_overhangs:
        best_right_overhang = select_best_overhang(
            contig_stats.right_overhangs, min_extension, max_homopolymer
        )

    # If no valid overhangs found, return None
    if not best_left_overhang and not best_right_overhang:
        return None

    # Apply extensions in order: left first, then right
    # This maintains correct positioning since left extensions don't affect right positions
    working_sequence = original_sequence
    final_extension_info = {}

    # Process left extension first
    if best_left_overhang:
        if not dry_run:
            try:
                working_sequence, left_ext_info = apply_contig_extension(
                    working_sequence, best_left_overhang, contig_stats.contig_length
                )
                # Update the extension info with left extension details
                final_extension_info.update(
                    {'left_' + k: v for k, v in left_ext_info.items()}
                )
                final_extension_info['has_left_extension'] = True
            except ValueError as e:
                warnings.append(f'Left extension failed for {contig_name}: {e}')
                best_left_overhang = None
        else:
            # Simulate left extension for dry run
            final_extension_info.update(
                {
                    'left_overhang_length': best_left_overhang.length,
                    'left_read_name': best_left_overhang.read_name,
                    'left_trim_length': 0,
                    'has_left_extension': True,
                }
            )

    # Process right extension second
    if best_right_overhang:
        if not dry_run:
            try:
                # For right extension after left extension, we need to adjust the overhang position
                # The right overhang alignment coordinates are relative to the original contig,
                # but we need to apply it to the extended sequence
                adjusted_right_overhang = OverhangInfo(
                    sequence=best_right_overhang.sequence,
                    length=best_right_overhang.length,
                    # Adjust alignment positions to account for left extension
                    alignment_pos=best_right_overhang.alignment_pos
                    + (len(working_sequence) - contig_stats.contig_length),
                    alignment_end=best_right_overhang.alignment_end
                    + (len(working_sequence) - contig_stats.contig_length),
                    read_name=best_right_overhang.read_name,
                    is_left=best_right_overhang.is_left,
                    clip_length=best_right_overhang.clip_length,
                    anchor_length=best_right_overhang.anchor_length,
                    contig_name=best_right_overhang.contig_name,
                )

                # Apply to the current working sequence (which may already include left extension)
                working_sequence, right_ext_info = apply_contig_extension(
                    working_sequence,
                    adjusted_right_overhang,
                    # Use the current length, not the original length
                    len(working_sequence),
                )
                # Update the extension info with right extension details
                final_extension_info.update(
                    {'right_' + k: v for k, v in right_ext_info.items()}
                )
                final_extension_info['has_right_extension'] = True
            except ValueError as e:
                warnings.append(f'Right extension failed for {contig_name}: {e}')
                best_right_overhang = None
        else:
            # Simulate right extension for dry run
            final_extension_info.update(
                {
                    'right_overhang_length': best_right_overhang.length,
                    'right_read_name': best_right_overhang.read_name,
                    'right_trim_length': 0,
                    'has_right_extension': True,
                }
            )

    # If all extensions failed, return None
    if not final_extension_info:
        return None

    # Add overall extension info
    final_extension_info.update(
        {
            'original_length': contig_stats.contig_length,
            'final_length': len(working_sequence)
            if not dry_run
            else contig_stats.contig_length
            + (best_left_overhang.length if best_left_overhang else 0)
            + (best_right_overhang.length if best_right_overhang else 0),
            'contig_name': contig_name,
        }
    )

    # For backward compatibility, set primary extension info to the first successful extension
    if best_left_overhang:
        final_extension_info.update(
            {
                'overhang_length': best_left_overhang.length,
                'read_name': best_left_overhang.read_name,
                'is_left': True,
                'trim_length': final_extension_info.get('left_trim_length', 0),
            }
        )
    elif best_right_overhang:
        final_extension_info.update(
            {
                'overhang_length': best_right_overhang.length,
                'read_name': best_right_overhang.read_name,
                'is_left': False,
                'trim_length': final_extension_info.get('right_trim_length', 0),
            }
        )

    # Count motifs if patterns provided
    motif_counts = {}
    if motif_patterns:
        target_seq = working_sequence if not dry_run else original_sequence

        # Add extensions to target sequence for dry run
        if dry_run:
            if best_left_overhang:
                target_seq = best_left_overhang.sequence + target_seq
            if best_right_overhang:
                target_seq = target_seq + best_right_overhang.sequence

        for pattern_name, pattern_str in motif_patterns.items():
            matches = re.findall(pattern_str, target_seq)
            motif_counts[pattern_name] = len(matches)

    return ExtensionResult(
        contig_name=contig_name,
        original_length=contig_stats.contig_length,
        extended_sequence=working_sequence,
        extension_info=final_extension_info,
        warnings=warnings,
        motif_counts=motif_counts,
    )
