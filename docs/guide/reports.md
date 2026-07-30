# Reading the report

`teloclip extend` can emit two reports, and they answer different questions.

| Option | Format | Answers |
|---|---|---|
| `--stats-report` | Markdown | *What changed?* |
| `--html-report` | Self-contained HTML | *Should I believe it?* |

The Markdown report is described first; the [HTML report](#the-html-report) is
covered at the end.

Pass `--stats-report` a path to write a file, `-` for stdout, or omit the option
and it goes to stderr interleaved with the log.

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

---

## The HTML report

`--html-report PATH` writes a single self-contained file: no CDN, no separate
stylesheet, nothing fetched at open time. It can be emailed, copied off a
cluster, or archived and still render years later.

It repeats the summary and extensions table, then adds two things the Markdown
report cannot express.

### Overhang depth across contigs

A strip plot with one point per contig end, blue for left and orange for right,
against a median reference line. Hover any point for the contig name and depth.

The shape is the message: a healthy assembly is a flat band near the median.
Points far above it are ringed, labelled with the contig name, and listed in the
anomalous-coverage section. A table view below carries the same numbers, so
nothing is reachable only by hovering.

### Overhang alignments

One scrollable block per contig end, showing every supporting read laid out
against the terminus:

```
                              │ terminus
contig                        │GGACCCTAACCCTAACCCTGCTAGATT…
▸ read_1  CCCTAACCCTAACCCTAACC│TAACCCTAACCCTAACCCTGCTAGATT…
  read_2  CCCTAACCCTAACCCTAACC│TAACCCTAACCCTAACCCTGCTAGATT…
  read_3        CCCTAACCCTAACC│TAACCCTAACCCTAACCCTGCTAGATT…
```

- **Grey** is the anchored portion of the read — the part that aligned.
- **Colour** is the soft clip, in that end's series colour.
- The **vertical rule** is the contig/overhang boundary.
- **`▸`** marks the read the extension actually used.
- Highlighted bases are motif matches, when `--count-motifs` was given.

Hover a row for the read name, how many bases it adds and how many it trims.

This view answers, in one place, the questions that otherwise require exporting
reads and opening an alignment viewer:

- **Do the clips agree?** Reads telling the same story is the evidence that the
  overhang is real. Reads disagreeing past the rule suggests non-specific
  alignment.
- **Does the anchor match the contig?** Compare the grey region against the
  `contig` row directly above it. Disagreement there means the read is not
  placed where it claims.
- **Is the repeat in the right place?** Motif highlighting should accumulate
  *past* the rule, not inside the contig.
- **Was the right read chosen?** The marked read is shown among the candidates
  it beat, with its trim cost visible.

Blocks open scrolled to the terminus. `--html-max-reads` (default 25) caps the
rows per end; the reads contributing the most sequence come first and the
remainder are counted in a note beneath. Clips are drawn to 300 bases, which is
well past the point of reading them base by base.
