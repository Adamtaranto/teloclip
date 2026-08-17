"""
Memory-efficient analysis functions for processing large genomes.

This module provides streaming analysis functions that process contigs
individually to minimize memory usage.
"""

from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional

import pysam

from .analysis import ContigStats, OverhangInfo
from .overhang import classify, ends_from_aligned_segment, overhang_info_from_call


def collect_contig_overhangs_streaming(
    bam_file: pysam.AlignmentFile,
    contig_name: str,
    contig_length: int,
    max_break: int = 10,
    min_clip: int = 1,
    min_anchor: int = 500,
    anchor_context: int = 0,
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
        Maximum tolerated gap between a contig terminus and the start of the
        alignment, inclusive (default: 10).
    min_clip : int, optional
        Minimum number of clipped bases required to lie past the terminus
        (default: 1).
    min_anchor : int, optional
        Minimum number of anchoring (M/=/X) bases required (default: 500).
    anchor_context : int, optional
        Aligned bases adjacent to each clip to retain for display
        (default: 0). Only the HTML report needs these.

    Returns
    -------
    ContigStats
        Statistics for this contig's overhangs.

    Notes
    -----
    Acceptance is delegated to :func:`teloclip.core.overhang.classify`, which applies
    the same rules used by ``filter`` and ``extract``. See that module for the
    coordinate convention.
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

        # Require a soft clip, and reject hard-clipped reads: their clipped
        # bases are absent from the record, so there is no sequence to graft.
        cigar_string = alignment.cigarstring
        if not cigar_string or 'S' not in cigar_string or 'H' in cigar_string:
            continue

        ends = ends_from_aligned_segment(alignment, contig_length, contig_name)
        left_call, right_call = classify(
            ends, max_break=max_break, min_clip=min_clip, min_anchor=min_anchor
        )

        # A read accepted at both ends spans the whole contig and hangs off
        # each side. That is worth surfacing: it usually means a very short
        # contig, or a circular molecule whose ends are the same locus.
        spans_both = left_call.accepted and right_call.accepted

        if left_call.accepted:
            info = overhang_info_from_call(ends, left_call, anchor_context)
            info.spans_both_ends = spans_both
            contig_stats.left_overhangs.append(info)

        if right_call.accepted:
            info = overhang_info_from_call(ends, right_call, anchor_context)
            info.spans_both_ends = spans_both
            contig_stats.right_overhangs.append(info)

    return contig_stats


def stream_contigs_for_extension(
    bam_file: pysam.AlignmentFile,
    contig_dict: Dict[str, int],
    min_overhangs: int = 1,
    max_break: int = 10,
    min_clip: int = 1,
    min_anchor: int = 500,
    anchor_context: int = 0,
) -> Iterator[tuple]:
    """
    Stream contigs that meet criteria for extension.

    Contigs are processed one at a time and never held collectively, so peak
    memory is set by the largest contig rather than by the assembly.

    Parameters
    ----------
    bam_file : pysam.AlignmentFile
        Opened BAM file with index.
    contig_dict : Dict[str, int]
        Dictionary mapping contig names to lengths.
    min_overhangs : int, optional
        Minimum number of overhangs required (default: 1).
    max_break : int, optional
        Maximum tolerated gap between a contig terminus and the alignment
        (default: 10).
    min_clip : int, optional
        Minimum number of clipped bases required past the terminus (default: 1).
    min_anchor : int, optional
        Minimum number of anchoring (M/=/X) bases required (default: 500).
    anchor_context : int, optional
        Aligned bases adjacent to each clip to retain for display (default: 0).

    Yields
    ------
    tuple
        (contig_name, contig_stats) for contigs that meet extension criteria.
    """
    for contig_name, contig_length in contig_dict.items():
        stats = collect_contig_overhangs_streaming(
            bam_file,
            contig_name,
            contig_length,
            max_break,
            min_clip,
            min_anchor,
            anchor_context,
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
    end_motif_counts: Dict[str, Dict[str, int]] = field(default_factory=dict)


def process_single_contig_extension(
    contig_name: str,
    contig_stats: ContigStats,
    original_sequence: str,
    min_extension: int = 1,
    max_homopolymer: int = 50,
    motif_patterns: Optional[Dict[str, str]] = None,
    dry_run: bool = False,
    terminal_length: int = 0,
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
    terminal_length : int, optional
        Number of original terminal bases to include in the per-end motif
        counting window, alongside the length of the extension at that end
        (default: 0).

    Returns
    -------
    Optional[ExtensionResult]
        Extension result if successful, None otherwise.
    """
    # Import here to avoid circular imports
    import re

    from .analysis import select_best_overhang
    from .extension import apply_contig_extension, calculate_extension_position

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
            # Simulate left extension for dry run. The trim must be reported
            # here too, otherwise a dry run predicts a different final length
            # than the real run produces for any contig needing trimming.
            _, left_trim = calculate_extension_position(
                best_left_overhang.alignment_pos,
                best_left_overhang.alignment_end,
                contig_stats.contig_length,
                True,
            )
            final_extension_info.update(
                {
                    'left_overhang_length': best_left_overhang.length,
                    'left_read_name': best_left_overhang.read_name,
                    'left_trim_length': left_trim,
                    'left_net_gain': best_left_overhang.length - left_trim,
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
                    net_gain=best_right_overhang.net_gain,
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
            # Simulate right extension for dry run (see the left branch above).
            _, right_trim = calculate_extension_position(
                best_right_overhang.alignment_pos,
                best_right_overhang.alignment_end,
                contig_stats.contig_length,
                False,
            )
            final_extension_info.update(
                {
                    'right_overhang_length': best_right_overhang.length,
                    'right_read_name': best_right_overhang.read_name,
                    'right_trim_length': right_trim,
                    'right_net_gain': best_right_overhang.length - right_trim,
                    'has_right_extension': True,
                }
            )

    # If all extensions failed, return None
    if not final_extension_info:
        return None

    # Add overall extension info. The dry-run length must account for trimming
    # as well as addition, so that a dry run and a real run agree.
    if dry_run:
        predicted_length = (
            contig_stats.contig_length
            + final_extension_info.get('left_net_gain', 0)
            + final_extension_info.get('right_net_gain', 0)
        )
    else:
        predicted_length = len(working_sequence)

    final_extension_info.update(
        {
            'original_length': contig_stats.contig_length,
            'final_length': predicted_length,
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
    end_motif_counts: Dict[str, Dict[str, int]] = {}
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

        # Count motifs in an explicit window at each end of the extended
        # sequence: the original terminal screening window plus whatever was
        # added at that end. This makes the counts directly comparable with the
        # pre-extension terminal counts. The window uses the net gain, since
        # that is how much longer the sequence actually is at that end.
        left_added = final_extension_info.get('left_net_gain', 0)
        right_added = final_extension_info.get('right_net_gain', 0)
        left_window = min(terminal_length + left_added, len(target_seq))
        right_window = min(terminal_length + right_added, len(target_seq))

        end_motif_counts = {'left': {}, 'right': {}}
        for pattern_name, pattern_str in motif_patterns.items():
            left_seq = target_seq[:left_window] if left_window > 0 else ''
            right_seq = target_seq[-right_window:] if right_window > 0 else ''
            end_motif_counts['left'][pattern_name] = len(
                re.findall(pattern_str, left_seq)
            )
            end_motif_counts['right'][pattern_name] = len(
                re.findall(pattern_str, right_seq)
            )

    return ExtensionResult(
        contig_name=contig_name,
        original_length=contig_stats.contig_length,
        extended_sequence=working_sequence,
        extension_info=final_extension_info,
        warnings=warnings,
        motif_counts=motif_counts,
        end_motif_counts=end_motif_counts,
    )
