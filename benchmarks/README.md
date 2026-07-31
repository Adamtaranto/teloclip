# Benchmarks

Performance benchmarks for teloclip's hot paths, run under
[CodSpeed](https://codspeed.io) via `pytest-codspeed`.

These live outside `tests/` deliberately. The test run is configured with
coverage instrumentation in `[tool.pytest.ini_options]`, which would distort the
measurements, and `testpaths = ["tests"]` keeps them out of a normal `pytest`
invocation.

## Running

```bash
pip install '.[test]'

# Measure
pytest benchmarks --codspeed --no-cov

# Just check they still run correctly, without measuring
pytest benchmarks --no-cov
```

CI runs them on every push to `main` and every pull request
(`.github/workflows/codspeed.yml`), which requires a `CODSPEED_TOKEN` secret.

## What is covered

| Area | Benchmarks | Why |
|---|---|---|
| `overhang` | CIGAR parsing, span/anchor reduction, `classify` | Runs once per alignment in every sub-command; the hottest path in the codebase |
| `filter` / `extract` | `processSamlines`, `EnhancedStreamingSamFilter`, with and without motifs | Consume the entire SAM stream, so they dominate wall-clock on whole-genome runs |
| `extend` | Per-contig overhang collection, terminal motif screening, overhang ranking | The extend inner loop, plus the regex scanning that grows with `--screen-terminal-bases` |

Literal and fuzzy motif patterns are benchmarked separately: fuzzy patterns
carry bounded repetition and are markedly more expensive to match, so mixing
them would hide a regression in either.

## Input data

The SAM benchmarks reuse `tests/data/test.sam` but repeat its records to roughly
25,000 alignments, so the measurement reflects the algorithm rather than fixture
overhead. The extend benchmarks use the indexed `tests/data/test.bam` and
`test.fna` directly.

## Adding a benchmark

Mark the function with `@pytest.mark.benchmark` and keep the measured body to
the work being studied — build inputs outside it. End with an assertion that the
work actually happened, so an accidentally empty input is a failure rather than
a suspiciously fast result.
