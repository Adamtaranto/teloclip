"""
Benchmarks for the filter and extract streaming paths.

``processSamlines`` and ``EnhancedStreamingSamFilter`` both consume an entire
SAM stream, so they dominate wall-clock for whole-genome runs. Both write to
stdout, which is redirected here so the measurement reflects the filtering work
rather than terminal or pipe I/O.
"""

import contextlib
import os

import pytest

from teloclip.motifs import make_fuzzy_motif_regex, make_motif_regex
from teloclip.samops import EnhancedStreamingSamFilter, processSamlines


@contextlib.contextmanager
def _discard_stdout():
    """
    Redirect stdout to the null device for the duration of the block.

    ``processSamlines`` writes kept records straight to stdout. Capturing into
    an in-memory buffer would make the benchmark partly a measurement of string
    concatenation, so the output is discarded instead.

    Yields
    ------
    None
    """
    with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
            yield


@pytest.mark.benchmark
def test_process_samlines_no_motifs(sam_lines, contigs):
    """Benchmark the filter path with clip criteria only."""
    with _discard_stdout():
        processSamlines(
            iter(sam_lines),
            contigs,
            motif_list=None,
            max_break=50,
            min_clip=1,
            min_anchor=100,
        )


@pytest.mark.benchmark
def test_process_samlines_with_motifs(sam_lines, contigs):
    """Benchmark the filter path with motif matching in clipped regions.

    Motif matching adds a regex scan per surviving read, and is the usual
    configuration in practice.
    """
    motifs = [make_motif_regex('TTAGGG'), make_motif_regex('CCCTAA')]

    with _discard_stdout():
        processSamlines(
            iter(sam_lines),
            contigs,
            motif_list=motifs,
            max_break=50,
            min_clip=1,
            min_anchor=100,
        )


@pytest.mark.benchmark
def test_process_samlines_fuzzy_motifs(sam_lines, contigs):
    """Benchmark the filter path with fuzzy motif patterns.

    Fuzzy patterns carry bounded repetition, which is markedly more expensive
    to match than the literal form.
    """
    motifs = [make_fuzzy_motif_regex('TTAGGG'), make_fuzzy_motif_regex('CCCTAA')]

    with _discard_stdout():
        processSamlines(
            iter(sam_lines),
            contigs,
            motif_list=motifs,
            max_break=50,
            min_clip=1,
            min_anchor=100,
        )


@pytest.mark.benchmark
def test_enhanced_streaming_filter(sam_lines, contigs):
    """Benchmark the extract path, which yields per-overhang dictionaries."""
    stream = EnhancedStreamingSamFilter(
        iter(sam_lines),
        contigs,
        max_break=50,
        min_clip=1,
        min_anchor=100,
    )

    # Drain the generator; the work is in producing the records, not consuming.
    count = sum(1 for _ in stream)

    # Guard against silently benchmarking an empty stream.
    assert count > 0


@pytest.mark.benchmark
def test_enhanced_streaming_filter_with_motifs(sam_lines, contigs):
    """Benchmark the extract path with per-overhang motif counting."""
    stream = EnhancedStreamingSamFilter(
        iter(sam_lines),
        contigs,
        max_break=50,
        min_clip=1,
        min_anchor=100,
        motif_patterns={'TTAGGG': 'TTAGGG', 'CCCTAA': 'CCCTAA'},
    )

    count = sum(1 for _ in stream)

    assert count > 0
