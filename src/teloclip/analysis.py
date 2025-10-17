"""
Overhang analysis and statistics collection for contig extension.

This module provides functionality to analyze soft-clipped alignments and collect
statistics about overhanging sequences that can be used to extend draft contigs.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Iterator
import statistics
from .samops import splitCIGAR, checkClips


@dataclass
class OverhangInfo:
    """Information about a single overhanging sequence."""

    sequence: str
    length: int
    alignment_pos: int
    alignment_end: int
    read_name: str
    is_left: bool
    clip_length: int
    anchor_length: int


@dataclass
class ContigStats:
    """Statistics for overhangs at both ends of a contig."""

    contig_name: str
    contig_length: int
    left_overhangs: List[OverhangInfo] = field(default_factory=list)
    right_overhangs: List[OverhangInfo] = field(default_factory=list)

    @property
    def left_count(self) -> int:
        """Number of left overhangs."""
        return len(self.left_overhangs)

    @property
    def right_count(self) -> int:
        """Number of right overhangs."""
        return len(self.right_overhangs)

    @property
    def left_total_length(self) -> int:
        """Total length of all left overhangs."""
        return sum(oh.length for oh in self.left_overhangs)

    @property
    def right_total_length(self) -> int:
        """Total length of all right overhangs."""
        return sum(oh.length for oh in self.right_overhangs)


def collect_overhang_stats(
    sam_lines: Iterator[str], contig_dict: Dict[str, int]
) -> Dict[str, ContigStats]:
    """
    Collect overhang statistics from SAM file lines.

    Parameters
    ----------
    sam_lines : Iterator[str]
        Iterator of SAM file lines
    contig_dict : Dict[str, int]
        Dictionary mapping contig names to their lengths

    Returns
    -------
    Dict[str, ContigStats]
        Dictionary mapping contig names to their overhang statistics
    """
    stats = {}

    # Initialize stats for all contigs
    for contig_name, contig_length in contig_dict.items():
        stats[contig_name] = ContigStats(contig_name, contig_length)

    for line in sam_lines:
        line = line.strip()
        if not line or line.startswith('@'):
            continue

        fields = line.split('\t')
        if len(fields) < 11:
            continue

        read_name = fields[0]
        flag = int(fields[1])
        ref_name = fields[2]
        pos = int(fields[3])
        cigar = fields[5]
        seq = fields[9]

        # Skip unmapped reads
        if flag & 4 or ref_name == '*' or ref_name not in contig_dict:
            continue

        # Skip secondary/supplementary alignments
        if flag & 0x900:
            continue

        # Analyze CIGAR for clips
        cigar_ops = splitCIGAR(cigar)
        left_clip, right_clip = checkClips(cigar)

        # Convert None to 0 for easier handling
        left_clip = left_clip or 0
        right_clip = right_clip or 0

        if left_clip == 0 and right_clip == 0:
            continue

        contig_length = contig_dict[ref_name]
        alignment_end = (
            pos + sum(length for length, op in cigar_ops if op in 'MDN=X') - 1
        )

        # Check for left clip at contig start
        if left_clip > 0 and pos <= 10:  # Allow some tolerance
            overhang_seq = seq[:left_clip]
            anchor_length = len(seq) - left_clip

            overhang = OverhangInfo(
                sequence=overhang_seq,
                length=left_clip,
                alignment_pos=pos,
                alignment_end=alignment_end,
                read_name=read_name,
                is_left=True,
                clip_length=left_clip,
                anchor_length=anchor_length,
            )
            stats[ref_name].left_overhangs.append(overhang)

        # Check for right clip at contig end
        if (
            right_clip > 0 and alignment_end >= contig_length - 10
        ):  # Allow some tolerance
            overhang_seq = seq[-right_clip:]
            anchor_length = len(seq) - right_clip

            overhang = OverhangInfo(
                sequence=overhang_seq,
                length=right_clip,
                alignment_pos=pos,
                alignment_end=alignment_end,
                read_name=read_name,
                is_left=False,
                clip_length=right_clip,
                anchor_length=anchor_length,
            )
            stats[ref_name].right_overhangs.append(overhang)

    return stats


def calculate_overhang_statistics(
    stats_dict: Dict[str, ContigStats],
) -> Dict[str, Dict[str, float]]:
    """
    Calculate statistical measures for overhang lengths across all contigs.

    Parameters
    ----------
    stats_dict : Dict[str, ContigStats]
        Dictionary of contig statistics

    Returns
    -------
    Dict[str, Dict[str, float]]
        Dictionary with 'left' and 'right' keys containing statistical measures
    """
    left_lengths = []
    right_lengths = []

    for contig_stats in stats_dict.values():
        left_lengths.extend([oh.length for oh in contig_stats.left_overhangs])
        right_lengths.extend([oh.length for oh in contig_stats.right_overhangs])

    def calc_stats(lengths: List[int]) -> Dict[str, float]:
        if not lengths:
            return {'mean': 0.0, 'median': 0.0, 'std_dev': 0.0, 'min': 0.0, 'max': 0.0}

        return {
            'mean': statistics.mean(lengths),
            'median': statistics.median(lengths),
            'std_dev': statistics.stdev(lengths) if len(lengths) > 1 else 0.0,
            'min': min(lengths),
            'max': max(lengths),
        }

    return {
        'left': calc_stats(left_lengths),
        'right': calc_stats(right_lengths),
        'combined': calc_stats(left_lengths + right_lengths),
    }


def identify_outlier_contigs(
    stats_dict: Dict[str, ContigStats], threshold: float = 2.0
) -> Dict[str, List[str]]:
    """
    Identify contigs with outlier overhang patterns using Z-score analysis.

    Parameters
    ----------
    stats_dict : Dict[str, ContigStats]
        Dictionary of contig statistics
    threshold : float, optional
        Z-score threshold for outlier detection (default: 2.0)

    Returns
    -------
    Dict[str, List[str]]
        Dictionary with 'left_outliers' and 'right_outliers' keys containing contig names
    """

    def calculate_z_scores(values: List[float]) -> List[float]:
        if len(values) <= 1:
            return [0.0] * len(values)

        mean_val = statistics.mean(values)

        try:
            std_val = statistics.stdev(values)
        except statistics.StatisticsError:
            return [0.0] * len(values)

        if std_val == 0:
            return [0.0] * len(values)

        return [(val - mean_val) / std_val for val in values]

    # Collect overhang counts per contig
    left_counts = []
    right_counts = []
    contig_names = []

    for contig_name, contig_stats in stats_dict.items():
        contig_names.append(contig_name)
        left_counts.append(contig_stats.left_count)
        right_counts.append(contig_stats.right_count)

    # Calculate Z-scores
    left_z_scores = calculate_z_scores(left_counts)
    right_z_scores = calculate_z_scores(right_counts)

    # Identify outliers
    left_outliers = []
    right_outliers = []

    for i, contig_name in enumerate(contig_names):
        if abs(left_z_scores[i]) > threshold:
            left_outliers.append(contig_name)
        if abs(right_z_scores[i]) > threshold:
            right_outliers.append(contig_name)

    return {'left_outliers': left_outliers, 'right_outliers': right_outliers}


def rank_overhangs_by_length(overhangs_list: List[OverhangInfo]) -> List[OverhangInfo]:
    """
    Sort overhangs by length in descending order.

    Parameters
    ----------
    overhangs_list : List[OverhangInfo]
        List of overhang information objects

    Returns
    -------
    List[OverhangInfo]
        Sorted list of overhangs (longest first)
    """
    return sorted(overhangs_list, key=lambda oh: oh.length, reverse=True)


def detect_homopolymer_runs(
    sequence: str, min_length: int = 50
) -> List[Tuple[str, int, int, int]]:
    """
    Detect homopolymer runs in a sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence to analyze
    min_length : int, optional
        Minimum length of homopolymer run to report (default: 50)

    Returns
    -------
    List[Tuple[str, int, int, int]]
        List of tuples (nucleotide, start_pos, end_pos, length) for each run
    """
    runs = []

    if not sequence:
        return runs

    current_base = sequence[0].upper()
    run_start = 0
    run_length = 1

    for i in range(1, len(sequence)):
        base = sequence[i].upper()

        if base == current_base:
            run_length += 1
        else:
            if run_length >= min_length:
                runs.append(
                    (current_base, run_start, run_start + run_length - 1, run_length)
                )

            current_base = base
            run_start = i
            run_length = 1

    # Check final run
    if run_length >= min_length:
        runs.append((current_base, run_start, run_start + run_length - 1, run_length))

    return runs


def select_best_overhang(
    overhangs: List[OverhangInfo], min_extension: int = 1, max_homopolymer: int = 50
) -> Optional[OverhangInfo]:
    """
    Select the best overhang for extension based on length and quality.

    Parameters
    ----------
    overhangs : List[OverhangInfo]
        List of available overhangs
    min_extension : int, optional
        Minimum overhang length required (default: 1)
    max_homopolymer : int, optional
        Maximum allowed homopolymer run length (default: 50)

    Returns
    -------
    Optional[OverhangInfo]
        Best overhang for extension, or None if no suitable overhang found
    """
    if not overhangs:
        return None

    # Filter by minimum length
    candidates = [oh for oh in overhangs if oh.length >= min_extension]

    if not candidates:
        return None

    # Sort by length (longest first)
    candidates = rank_overhangs_by_length(candidates)

    # Check for homopolymer runs in order of preference
    for overhang in candidates:
        homo_runs = detect_homopolymer_runs(overhang.sequence, max_homopolymer)
        if not homo_runs:  # No concerning homopolymer runs
            return overhang

    # If all have homopolymer issues, return the longest one anyway
    # but this should be flagged in reporting
    return candidates[0] if candidates else None
