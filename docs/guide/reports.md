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

It repeats the summary and extensions table, then adds things the Markdown
report cannot express: three charts and the read-level alignments.

All three charts share a convention. Blue is the left end of a contig, orange
the right. Marks are keyboard-focusable and show the same detail on focus as on
hover. **Clicking any mark selects that contig in all three charts at once**, so
an end that looks unremarkable in one view can be found in the others; the rest
recede rather than disappearing, because the distribution is the context for the
selection. Press <kbd>Escape</kbd>, or click the background of any chart, to
clear. Every chart has a table view below it, so no value is reachable only by
hovering.

### Overhang depth across contigs

A strip plot with one point per contig end against a median reference line.

The shape is the message: a healthy assembly is a flat band near the median.
Points far above it are ringed, labelled with the contig name, and listed in the
anomalous-coverage section.

### Overhang length distribution

Depth says how many reads support an end. This says how far they reach past it,
which is what determines how much sequence an extension can actually recover — a
tight distribution at a few hundred bases reads very differently from a broad one
running into kilobases, and neither is visible in a count.

Each contig is one shape, split at its centre line: the left half is the left end,
the right half the right end. Width is the proportion of reads at that length, the
dark tick is the median, and the short rule beside it is the interquartile range.

Two things are worth knowing about how it is drawn:

- **Ends with fewer than four reads are drawn as individual points**, not as a
  curve. A distribution estimated from three reads describes the estimator rather
  than the data.
- **The length axis is clipped at the 99th percentile** so that a single read
  running far past a contig end cannot flatten every other shape into a line. When
  that happens the caption says so and names the true longest overhang, so the
  figure you are looking for is never hidden by the clip.

### Overhang depth against length

One point per contig end, depth against median overhang length. This is where the
two preceding charts are reconciled, and where the corner a point lands in tells
you what kind of end it is:

| Position | Reading |
|---|---|
| Top right — deep and long | A collapsed repeat or rDNA array at the terminus. Extending from a single read here is rarely meaningful. |
| Top left — deep but short | Reads drawn in from elsewhere, as an organellar contig or a collapsed repeat does. |
| Bottom right — long but shallow | The ordinary shape of a telomere the assembly stopped short of. This is the case extension exists to serve. |
| Bottom left — shallow and short | Little evidence either way. Not an anomaly, just a quiet end. |

Hovering a point gives the full figures for that end: contig length, depth,
median and longest overhang, the best net gain available, and which anomaly
flags apply.

### Overhang alignments

One scrollable block per contig end, showing every supporting read laid out
against the terminus:

```text
          -20        0        +20
           |         |         |
contig     ACGTACGT--ACCCTAACCCTGCTAGATT…
▸ read_1  CCCTAACCCT--ACCCTAACCCTGCTAGATT…
  read_2  CCCTAACCCTNNACCCTAACCCTGCTAGATT…
  read_3      TAACCCT--ACCC--ACCCTGCTAGATT…
```

- **Grey** is the anchored portion of the read — the part that aligned.
- **Colour** is the soft clip, in that end's series colour.
- The **vertical rule** is the contig/overhang boundary.
- The **ruler** across the top is in contig bases with `0` at the terminus,
  negative outside the contig and positive inside it.
- **`▸`** marks the read the extension actually used.
- **`⇆`** marks a read that overhangs *both* ends of the contig, which usually
  means a very short contig or a circular molecule.
- Highlighted bases are motif matches, when `--count-motifs` was given.

Reads are placed by walking their CIGAR, so an indel shifts a read against the
reference exactly as the aligner recorded it. A deletion leaves the read gapped
(`-`) where the contig still has bases; an insertion opens a column in *every*
other row, including the contig, so the rows below stay aligned. Pasting read
sequence beside the reference instead would drift by one column per indel and
make a good alignment look mismatched.

Hovering a row gives a table of the SAM record — alignment span, CIGAR, FLAG,
MAPQ, anchor and clip lengths — alongside what the read contributes and what it
costs in trimming.

The anchored region shown is at least 500 bases, or the longest anchor among
the reads at that end if that is longer, capped at 2,000 and never more than the
contig itself.

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

### Provenance

The report closes with the teloclip version, the time it was generated, and the
full command that produced it. A report is an artefact people keep and share, so
it has to be traceable back to the run behind it without relying on where the
file happens to sit.
