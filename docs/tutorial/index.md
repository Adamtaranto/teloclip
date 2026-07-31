# Tutorial

This tutorial walks through a complete teloclip run, from installing the tool to
validating the extended assembly.

Each stage is a separate page and builds on the previous one. The commands use a
consistent set of filenames so you can follow along end to end.

| Stage | Page | Produces |
|---|---|---|
| 1 | [Installation](install.md) | A working `teloclip` and `samtools` |
| 2 | [Mapping reads](mapping.md) | `raw.bam` and its index |
| 3 | [Filtering overhangs](filter.md) | `overhangs.bam` |
| 4 | [Extracting reads](extract.md) | Per-contig-end FASTA files |
| 5 | [Extending contigs](extend.md) | `extended.fasta` and a stats report |
| 6 | [Validating the result](validate.md) | Confidence, or a reason to stop |

## What you need

- A draft assembly in FASTA format (`ref.fa`)
- The long reads used to build it (`reads.fq.gz`)
- `minimap2` and `samtools`

Reads should be the raw long reads, not the corrected or trimmed set — teloclip
needs the sequence the assembler discarded, which correction may remove.

## Files used throughout

```
ref.fa              draft assembly
ref.fa.fai          its index
reads.fq.gz         raw long reads
raw.bam             all alignments
overhangs.bam       alignments clipped at contig ends
extended.fasta      the result
extension_report.md what changed and why
```

## Before you start

Know your organism's telomeric repeat. The canonical vertebrate repeat is
`TTAGGG`, but plants commonly use `TTTAGGG`, and many fungi and protists differ
again. Passing the wrong motif to `--motifs` will silently discard the reads you
want.

If you do not know it, run [`filter`](filter.md) without `--motifs` first and
inspect the clipped sequences with [`extract`](extract.md) — the repeat is
usually obvious once you look at them.
