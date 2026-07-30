"""
Overhang analysis and statistics collection for contig extension.

This module provides functionality to analyze soft-clipped alignments and collect
statistics about overhanging sequences that can be used to extend draft contigs.
"""

from dataclasses import dataclass, field
import logging
import statistics
from typing import Dict, List, Optional, Tuple


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
    contig_name: str
    # Bases the contig grows by if this overhang is applied: the clipped bases
    # lying past the terminus, i.e. clip length minus the gap that gets trimmed.
    # See teloclip.overhang for the full coordinate convention.
    net_gain: int = 0


@dataclass
class ContigStats:
    """Statistics for overhangs at both ends of a contig."""

    contig_name: str
    contig_length: int
    left_overhangs: List[OverhangInfo] = field(default_factory=list)
    right_overhangs: List[OverhangInfo] = field(default_factory=list)

    @property
    def left_count(self) -> int:
        """
        Number of left overhangs.

        Returns
        -------
        int
            The number of left overhangs (i.e., the length of self.left_overhangs).
        """
        return len(self.left_overhangs)

    @property
    def right_count(self) -> int:
        """
        Number of right overhangs.

        Returns
        -------
        int
            The number of right overhangs (i.e., the length of self.right_overhangs).
        """
        return len(self.right_overhangs)

    @property
    def left_total_length(self) -> int:
        """
        Total length of all left overhangs.

        Compute the sum of the lengths of every overhang in the
        ``left_overhangs`` sequence. Each element in that sequence is expected
        to expose a numeric ``length`` attribute. If the sequence is empty,
        the result is 0.

        Returns
        -------
        int
            Sum of the ``length`` attributes of all left overhangs.
        """
        return sum(oh.length for oh in self.left_overhangs)

    @property
    def right_total_length(self) -> int:
        """
        Total length of all right overhangs.

        Compute the sum of the ``length`` attribute for each overhang in
        ``self.right_overhangs``.

        Returns
        -------
        int
            Total length of all right overhangs. Returns 0 when ``self.right_overhangs`` is empty.

        Notes
        -----
        Each element of ``self.right_overhangs`` is expected to expose a numeric
        ``length`` attribute (typically an integer). This property does not
        modify the underlying collection.
        """
        return sum(oh.length for oh in self.right_overhangs)


def calculate_overhang_statistics(
    stats_dict: Dict[str, ContigStats],
) -> Dict[str, Dict[str, float]]:
    """
    Calculate statistical measures for overhang lengths across all contigs.

    Parameters
    ----------
    stats_dict : Dict[str, ContigStats]
        Dictionary of contig statistics.

    Returns
    -------
    Dict[str, Dict[str, float]]
        Dictionary with 'left' and 'right' keys containing statistical measures.
    """
    left_lengths = []
    right_lengths = []

    for contig_stats in stats_dict.values():
        left_lengths.extend([oh.length for oh in contig_stats.left_overhangs])
        right_lengths.extend([oh.length for oh in contig_stats.right_overhangs])

    def calc_stats(lengths: List[int]) -> Dict[str, float]:
        """
        Compute basic descriptive statistics for a sequence of lengths.

        Parameters
        ----------
        lengths : list of int
            Sequence of lengths to summarize. If empty, all returned statistics are 0.0.

        Returns
        -------
        dict
            Dictionary mapping statistic names to their values. Contains the following keys:
            - 'mean' : float
                Arithmetic mean of the values (0.0 if lengths is empty).
            - 'median' : float
                Median value (0.0 if lengths is empty).
            - 'std_dev' : float
                Sample standard deviation (uses statistics.stdev). Returns 0.0 if fewer than two values.
            - 'min' : float
                Minimum value (0.0 if lengths is empty).
            - 'max' : float
                Maximum value (0.0 if lengths is empty).

        Notes
        -----
        The function guards against statistics.StatisticsError by returning 0.0 for the standard
        deviation when the input has fewer than two elements. For an empty input, all statistics
        are defined to be 0.0 for convenience.

        Examples
        --------
        >>> calc_stats([1, 2, 3])
        {'mean': 2.0, 'median': 2.0, 'std_dev': 1.0, 'min': 1, 'max': 3}

        >>> calc_stats([])
        {'mean': 0.0, 'median': 0.0, 'std_dev': 0.0, 'min': 0.0, 'max': 0.0}
        """
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


# Below this many contigs, the spread of overhang counts carries too little
# information to call anything anomalous, so flagging declines to judge.
MIN_CONTIGS_FOR_ANOMALY_FLAGGING = 8

# Scale factors making each absolute-deviation statistic a consistent estimator
# of the standard deviation for normally distributed data.
_MAD_TO_SIGMA = 1.4826
_MEAN_AD_TO_SIGMA = 1.253314


def _modified_z_scores(values: List[float]) -> List[float]:
    """
    Score values by their deviation from the median, in robust sigma units.

    Uses the median and the median absolute deviation rather than the mean and
    standard deviation. The mean and standard deviation are themselves dragged
    upward by the very outliers being looked for, which masks them; the median
    and MAD are not.

    The MAD collapses to zero whenever more than half the values are identical,
    which is common here: most assemblies have a large majority of contigs
    carrying the same small number of terminal overhangs. In that case the mean
    absolute deviation is used as the scale instead, which is less robust but
    still resistant enough, and which is what lets a genuinely extreme contig be
    seen against a flat background.

    Parameters
    ----------
    values : List[float]
        Values to score. Order is preserved in the output.

    Returns
    -------
    List[float]
        One score per input value. All zeros when there are fewer than two
        values, or when every value is identical and no spread exists at all.
    """
    if len(values) <= 1:
        return [0.0] * len(values)

    median = statistics.median(values)
    deviations = [abs(v - median) for v in values]

    scale = statistics.median(deviations) * _MAD_TO_SIGMA
    if scale == 0:
        scale = statistics.fmean(deviations) * _MEAN_AD_TO_SIGMA
    if scale == 0:
        # Every value is identical: no spread, nothing to flag.
        return [0.0] * len(values)

    return [(v - median) / scale for v in values]


def flag_anomalous_overhang_coverage(
    stats_dict: Dict[str, ContigStats],
    threshold: float = 3.5,
    min_contigs: int = MIN_CONTIGS_FOR_ANOMALY_FLAGGING,
) -> Dict[str, List[str]]:
    """
    Flag contigs whose overhang counts are unusually high.

    A contig accumulating far more clipped reads at an end than its peers is
    worth a look: it often marks a collapsed repeat, a rDNA array, or an
    organellar contig pulling in reads from across the assembly. Extending such
    a contig from a single read is rarely meaningful.

    Only the high tail is flagged. A contig with unusually *few* overhangs is
    simply one with little evidence, which is not an anomaly and never a reason
    to withhold it from extension.

    This reports; it does not exclude. Whether a flagged contig should be left
    out is a judgement about the assembly that belongs to the user, who can act
    on it with ``--exclude-contigs``.

    Parameters
    ----------
    stats_dict : Dict[str, ContigStats]
        Overhang statistics per contig.
    threshold : float, optional
        Modified z-score above which a contig is flagged (default: 3.5, the
        conventional cutoff for this statistic).
    min_contigs : int, optional
        Minimum number of contigs required before flagging is attempted
        (default: 8).

    Returns
    -------
    Dict[str, List[str]]
        Contig names under the keys ``left_outliers`` and ``right_outliers``.
        Both are empty when there are too few contigs to judge.
    """
    contig_names = list(stats_dict)

    if len(contig_names) < min_contigs:
        return {'left_outliers': [], 'right_outliers': []}

    left_counts = [stats_dict[name].left_count for name in contig_names]
    right_counts = [stats_dict[name].right_count for name in contig_names]

    left_scores = _modified_z_scores(left_counts)
    right_scores = _modified_z_scores(right_counts)

    return {
        'left_outliers': [
            name for name, score in zip(contig_names, left_scores) if score > threshold
        ],
        'right_outliers': [
            name for name, score in zip(contig_names, right_scores) if score > threshold
        ],
    }


# Candidates whose net gain falls within this margin of the best available are
# treated as equally good, and the one trimming least assembly wins. The margin
# is the larger of an absolute floor and a fraction of the best gain, so it stays
# meaningful for both short and long overhangs.
GAIN_MARGIN_BASES = 10
GAIN_MARGIN_FRACTION = 0.05


def rank_overhangs_by_gain(
    overhangs_list: List[OverhangInfo],
    margin_bases: int = GAIN_MARGIN_BASES,
    margin_fraction: float = GAIN_MARGIN_FRACTION,
) -> List[OverhangInfo]:
    """
    Sort overhangs by the sequence they would add, preferring less trimming.

    Ranking is on ``net_gain`` rather than raw clip length. Applying an overhang
    trims the contig bases the read did not cover before grafting on the clip,
    so a long clip anchored well inside the contig can contribute less novel
    sequence than a shorter clip flush with the terminus.

    Maximising net gain alone is not quite right either, because trimming
    discards polished assembly consensus and replaces it with a single raw
    read's version of the same region. Gaining one extra base at the cost of
    thirteen bases of consensus is a bad trade. Candidates are therefore grouped
    into bands of comparable net gain, and within a band the one that trims
    least wins. Only a materially larger gain, one that clears the margin,
    outranks a candidate that leaves the contig intact.

    Parameters
    ----------
    overhangs_list : List[OverhangInfo]
        List of overhang information objects.
    margin_bases : int, optional
        Absolute floor on the band width, in bases (default: 10).
    margin_fraction : float, optional
        Band width as a fraction of the best available net gain (default: 0.05).

    Returns
    -------
    List[OverhangInfo]
        Sorted list of overhangs, best candidate first.
    """
    if not overhangs_list:
        return []

    best_gain = max(oh.net_gain for oh in overhangs_list)
    # At least 1, so the band arithmetic below cannot divide by zero.
    margin = max(1, margin_bases, int(best_gain * margin_fraction))

    def sort_key(oh: OverhangInfo) -> Tuple[int, int, int, int]:
        """
        Build the ranking key for one overhang.

        Parameters
        ----------
        oh : OverhangInfo
            Overhang to rank.

        Returns
        -------
        Tuple[int, int, int, int]
            ``(band, trim, -net_gain, -length)``, ascending. Banding by distance
            from the best gain keeps the ordering a proper total order, rather
            than the intransitive mess a pairwise "within margin" test would
            produce.
        """
        band = (best_gain - oh.net_gain) // margin
        # The gap that gets trimmed is whatever part of the clip did not clear
        # the terminus.
        trim = oh.length - oh.net_gain
        return (band, trim, -oh.net_gain, -oh.length)

    return sorted(overhangs_list, key=sort_key)


def detect_homopolymer_runs(
    sequence: str, min_length: int = 50
) -> List[Tuple[str, int, int, int]]:
    """
    Detect homopolymer runs in a sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence to analyze.
    min_length : int, optional
        Minimum length of homopolymer run to report (default: 50).

    Returns
    -------
    List[Tuple[str, int, int, int]]
        List of tuples (nucleotide, start_pos, end_pos, length) for each run.
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
        List of available overhangs.
    min_extension : int, optional
        Minimum number of novel bases the overhang must contribute (default: 1).
    max_homopolymer : int, optional
        Maximum allowed homopolymer run length (default: 50).

    Returns
    -------
    Optional[OverhangInfo]
        Best overhang for extension, or None if no suitable overhang found.
    """
    if not overhangs:
        return None

    # Filter on the sequence actually contributed, not the raw clip length.
    candidates = [oh for oh in overhangs if oh.net_gain >= min_extension]

    if not candidates:
        return None

    # Sort by net gain (most novel sequence first)
    candidates = rank_overhangs_by_gain(candidates)

    # Check for homopolymer runs in order of preference
    excluded_count = 0
    for overhang in candidates:
        homo_runs = detect_homopolymer_runs(overhang.sequence, max_homopolymer)
        if homo_runs:  # Has concerning homopolymer runs
            side = 'left' if overhang.is_left else 'right'
            warning_msg = (
                f'Excluding overhang from read {overhang.read_name} on {side} side of '
                f'contig {overhang.contig_name}: homopolymer run {homo_runs[0][0]}x{homo_runs[0][3]} '
                f'at position {homo_runs[0][1]}-{homo_runs[0][2]} (length {homo_runs[0][3]}) '
                f'exceeds threshold ({max_homopolymer})'
            )
            logging.warning(warning_msg)
            excluded_count += 1
        else:
            return overhang

    # If all candidates have homopolymer issues, return None
    if excluded_count > 0:
        logging.warning(
            f'All {excluded_count} candidate overhangs contain homopolymer runs '
            f'exceeding threshold ({max_homopolymer}). No extension will be applied.'
        )

    return None
