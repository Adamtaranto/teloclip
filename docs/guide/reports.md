# Reading the report

`teloclip extend --stats-report` writes a Markdown report. It renders as a web
page but is also aligned plain text, so it reads fine in a terminal.

Pass a path to write a file, `-` for stdout, or omit the option and it goes to
stderr interleaved with the log.

## Summary

```
| Metric                        |             Value |
| Mode                          | extension applied |
| Contigs in assembly           |                 1 |
| Contigs with overhang support |                 1 |
| Contigs extended              |                 1 |
| Contig ends extended          |            2 of 2 |
| Net bases gained              |               373 |
| Bases trimmed back            |                 2 |
| Raw overhang bases grafted    |               375 |
| Motifs counted                |    CCCTAA, TTAGGG |
| Total motif gain              |               +61 |
| Terminal screening window     |          1,000 bp |
| Contigs excluded              |                 0 |
| Warnings                      |                 0 |
```

The three length figures are distinct and worth reading together:

- **Net bases gained** — how much longer the assembly actually is
- **Bases trimmed back** — assembly consensus discarded to graft the clips on
- **Raw overhang bases grafted** — clipped read bases added, the sum of the two

**Contig ends extended** is the quickest sanity check. Compare it against twice
your chromosome count.

## Extensions Applied

```
| Contig | Original bp | Final bp | Left +bp | Right +bp | Total +bp | ... | Left trim | Right trim |
| chr1   |   6,030,967 | 6,031,340|     +160 |      +213 |      +373 | ... |         2 |          0 |
```

The `+bp` columns are **net**: the overhang grafted at that end, less the contig
bases trimmed to make room. So `Original bp + Total +bp = Final bp` on every
row — an identity you can check against the sequences.

The trim columns show what was given up. A large trim for a small gain deserves
a look; teloclip's ranking avoids that trade, but the numbers let you confirm.

The read name columns identify exactly which read supplied each extension, so
you can pull it from the BAM and inspect it.

## Telomere Motif Analysis

Only present with `--count-motifs`.

```
| Contig | End   | Motif  | Window bp |  Pre | Post | Gain |
| chr1   | left  | CCCTAA |     1,162 |    2 |   29 |  +27 |
| chr1   | right | TTAGGG |     1,213 |    1 |   35 |  +34 |
```

`Pre` counts the `--screen-terminal-bases` window of the original contig; `Post`
counts that same window plus what was added at that end, so the two are
comparable.

The window is reported as actually used. For contigs shorter than twice the
requested length it is halved, and the table shows the reduced figure.

Expect gain at the end where that motif belongs: `CCCTAA` at the 5' end and
`TTAGGG` at the 3' end for the canonical vertebrate repeat. Gain of zero is a
legitimate result — it means the extension added no telomeric sequence.

## Per-Contig Overhang Support

Read counts and the longest overhang seen at each end.

The longest overhang is not always the one used. A longer clip anchored further
inside the contig can contribute less novel sequence, and clips containing long
homopolymer runs are rejected. Where `Longest` exceeds the applied extension,
one of those applied.

## Contigs With Anomalous Overhang Coverage

Present only when something was flagged.

These ends carry far more clipped reads than the rest of the assembly, computed
with a median/MAD statistic on the high tail. That pattern usually means a
collapsed repeat, an rDNA array, or an organellar contig attracting reads from
across the assembly — cases where extending from one read is not meaningful.

**They are not excluded.** Whether extension is appropriate is a judgement about
your assembly. Review them, and re-run with `--exclude-contigs` if you agree.

Flagging is skipped on assemblies with fewer than 8 contigs, where the spread
carries too little information to call anything anomalous. It reports nothing
rather than guessing.

## Excluded Contigs

Contigs you named with `--exclude-contigs` or `--exclude-contigs-file`. They
still appear in the output FASTA, unmodified and in their original position.

## Warnings

Per-contig problems: extensions rejected for insufficient net gain, overhangs
rejected for homopolymer content, contigs in the index but missing from the
FASTA. Worth reading even when the run succeeds.
