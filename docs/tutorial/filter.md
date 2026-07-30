# 3. Filtering overhangs

`teloclip filter` reads SAM on stdin and writes SAM on stdout, keeping only
alignments that are soft-clipped at a contig end.

## Basic run

```bash
samtools view -h raw.bam \
  | teloclip filter --ref-idx ref.fa.fai \
  | samtools sort -o overhangs.bam
samtools index overhangs.bam
```

## Requiring telomeric motifs

Most overhangs at contig ends are not telomeres. Requiring the clipped region to
contain your organism's telomeric repeat removes most of the noise:

```bash
samtools view -h raw.bam \
  | teloclip filter \
      --ref-idx ref.fa.fai \
      --motifs TTAGGG \
      --min-repeats 3 \
  | samtools sort -o overhangs.bam
samtools index overhangs.bam
```

The reverse complement is searched automatically, so `--motifs TTAGGG` also
matches `CCCTAA` at the 5' end. Pass `--no-rev` to disable that.

`--fuzzy` tolerates ±1 base in homopolymer runs (`TTAGGG` becomes
`T{1,3}AG{2,4}`), which helps with ONT data where homopolymer length is
unreliable.

## The options that matter

| Option | Default | What it controls |
|---|---|---|
| `--max-break` | 50 | How far the alignment may stop short of the contig end |
| `--min-clip` | 1 | How far the clip must extend **past** the contig end |
| `--min-anchor` | 100 | Aligned (`M`/`=`/`X`) bases required, to trust the placement |
| `--min-repeats` | 1 | Consecutive motif copies required |
| `--match-anywhere` | off | Allow the motif anywhere in the read, not just the clip |

See [How overhangs are defined](../guide/overhangs.md) for exactly what
`--max-break` and `--min-clip` measure. In short: both are distances from the
contig terminus, `--max-break` inward and `--min-clip` outward.

!!! tip "Start permissive"

    Run without `--motifs` first and look at what you get. Tightening later is
    easy; discovering you filtered away the real telomeres is not.

## Reading the summary

`filter` reports what it discarded and why:

```
Processed 5482 SAM records.
Found 214 alignments soft-clipped at contig ends.
Found 3 alignments spanning entire contigs.
Output 187 alignments containing motif matches.
Exclusion summary:
  - Unmapped reads: 112
  - Secondary alignments: 0
  - No usable soft clip: 4903
  - Below min_anchor threshold (100bp): 41
  - Beyond max_break threshold (50bp): 168
  - Clip does not reach 1bp past contig end: 30
  - No telomeric motifs: 27
Total discarded: 5281 alignments after all filtering.
```

The buckets are mutually exclusive and sum to the discard total, so the figures
can be checked against each other.

It also logs per-contig counts and a histogram of overhang depth per contig end:

```
Overhang reads per contig end (24 ends):
    0-8  | ######################## 14
    9-17 | ########                  5
   18-26 | ###                       2
   ...
  90-100 | #                         1
```

A long right tail is worth investigating — see
[Collapsed repeats](../guide/biology.md#collapsed-repeats-and-rdna-arrays).

## Diagnosing an empty result

| Symptom | Likely cause |
|---|---|
| "No usable soft clip" is nearly everything | Normal — most reads are internal. Only worry if it is *literally* everything |
| Everything excluded by `max_break` | Alignments stop well short of contig ends. Raise `--max-break`, or your contigs have low-coverage tails |
| Everything excluded by `min_anchor` | Reads too short or too poorly aligned. Lower `--min-anchor` |
| Everything excluded by "No telomeric motifs" | Wrong motif for your organism, or try `--fuzzy` |
| Zero soft-clipped reads at all | Your aligner is hard-clipping. See [Mapping](mapping.md) |

Keep a copy of the log with `--logfile filter.log`.

Next: [Extracting reads](extract.md).
