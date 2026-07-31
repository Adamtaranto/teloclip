# Synthetic test data generators

These scripts build the synthetic dataset used by
`tests/integration/test_teloclip_extend_integration.py`. They are **not** tests:
they have no test functions and are not collected by pytest. They previously
lived under `tests/integration/` with `test_` prefixes, which meant pytest
collected three files contributing nothing.

The generated files are committed under `tests/integration/test_data/`, so you
only need these when regenerating or extending the dataset.

## Usage

Run from the repository root:

```bash
# 1. Synthetic contigs covering the edge cases the integration tests exercise
python scripts/testdata/generate_data.py

# 2. Alignments against those contigs, with overhangs at known positions
python scripts/testdata/generate_alignments.py

# 3. End-to-end check that the generated data behaves as the tests expect
python scripts/testdata/validate_integration.py [--quick] [--verbose]
```

`generate_alignments.py` writes SAM; convert and index it with samtools to
produce the `.bam`/`.bai` the integration tests load:

```bash
samtools sort tests/integration/test_data/synthetic_alignments.sam \
  -o tests/integration/test_data/synthetic_alignments_sorted.bam
samtools index tests/integration/test_data/synthetic_alignments_sorted.bam
```

## Note on expected values

Contig lengths are duplicated between `generate_data.py` and
`generate_alignments.py`; they must agree, or the alignments will place
overhangs at the wrong offsets. If you change one, change the other.
