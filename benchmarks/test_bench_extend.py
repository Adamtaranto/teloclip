"""
Benchmarks for the extend path.

Contig extension is dominated by two things: fetching and classifying the
alignments at each contig end, and scanning terminal windows for motifs. The
latter grows with ``--screen-terminal-bases`` and is easy to regress by
changing how the regex is applied.
"""

import pysam
import pytest

from teloclip.commands.extend import count_terminal_motifs, get_motif_regex
from teloclip.core.analysis import ContigStats, rank_overhangs_by_gain
from teloclip.core.overhang import (
    AlignmentEnds,
    classify_end,
    overhang_info_from_call,
)
from teloclip.core.streaming_analysis import collect_contig_overhangs_streaming


@pytest.mark.benchmark
def test_collect_contig_overhangs_streaming(contigs, bam_path):
    """Benchmark per-contig overhang collection from an indexed BAM.

    This is the extend command's inner loop: fetch each contig's alignments and
    classify both ends of every one.
    """
    with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        for _ in range(20):
            for name, length in contigs.items():
                collect_contig_overhangs_streaming(
                    bam, name, length, max_break=50, min_clip=1, min_anchor=100
                )


@pytest.mark.benchmark
def test_count_terminal_motifs(contigs, fasta_path):
    """Benchmark terminal-window motif screening across the assembly."""
    patterns = get_motif_regex('TTAGGG,CCCTAA', fuzzy=False)

    for _ in range(20):
        count_terminal_motifs(fasta_path, contigs, patterns, terminal_length=200)


@pytest.mark.benchmark
def test_count_terminal_motifs_fuzzy(contigs, fasta_path):
    """Benchmark terminal screening with fuzzy patterns.

    Bounded repetition makes these substantially more expensive than the
    literal patterns, so they are tracked separately.
    """
    patterns = get_motif_regex('TTAGGG,CCCTAA', fuzzy=True)

    for _ in range(20):
        count_terminal_motifs(fasta_path, contigs, patterns, terminal_length=200)


@pytest.mark.benchmark
def test_rank_overhangs_by_gain():
    """Benchmark overhang ranking at a realistic depth.

    A high-coverage contig end can accumulate hundreds of candidate overhangs,
    all of which are ranked before one is chosen.
    """
    contig_length = 100_000
    overhangs = []
    for i in range(500):
        ends = AlignmentEnds(
            contig_name='contig',
            contig_length=contig_length,
            aln_start=1 + (i % 40),
            aln_end=contig_length - (i % 7),
            left_clip=50 + (i % 300),
            right_clip=40 + (i % 200),
            anchor=1000,
            read_name=f'read{i}',
            sequence='ACGT' * 200,
        )
        call = classify_end(ends, True, max_break=50, min_clip=1)
        if call.accepted:
            overhangs.append(overhang_info_from_call(ends, call))

    for _ in range(200):
        rank_overhangs_by_gain(overhangs)


@pytest.mark.benchmark
def test_contig_stats_aggregation():
    """Benchmark the per-contig statistics properties.

    These are recomputed rather than cached, and are read repeatedly while
    building the report.
    """
    stats = ContigStats('contig', 100_000)
    ends = AlignmentEnds(
        contig_name='contig',
        contig_length=100_000,
        aln_start=1,
        aln_end=100_000,
        left_clip=120,
        right_clip=90,
        anchor=5000,
        read_name='read',
        sequence='ACGT' * 100,
    )
    left = classify_end(ends, True, max_break=50, min_clip=1)
    right = classify_end(ends, False, max_break=50, min_clip=1)
    stats.left_overhangs = [overhang_info_from_call(ends, left)] * 500
    stats.right_overhangs = [overhang_info_from_call(ends, right)] * 500

    total = 0
    for _ in range(2000):
        total += stats.left_count + stats.right_count
        total += stats.left_total_length + stats.right_total_length

    assert total > 0
