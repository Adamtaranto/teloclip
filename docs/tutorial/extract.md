# 4. Extracting reads

`teloclip extract` writes the overhang reads to one FASTA file per contig end,
so you can look at them before committing to an extension.

This step is optional but strongly recommended. It is the only stage where you
see the actual sequence teloclip proposes to add.

## Basic run

```bash
samtools view -h overhangs.bam \
  | teloclip extract \
      --ref-idx ref.fa.fai \
      --extract-dir split_by_contig
```

This writes files named `<contig>_L.fasta` and `<contig>_R.fasta` for the 5' and
3' ends respectively. The output directory is created if it does not exist.

You can pipe straight from `filter` without an intermediate file:

```bash
samtools view -h raw.bam \
  | teloclip filter --ref-idx ref.fa.fai --motifs TTAGGG \
  | teloclip extract --ref-idx ref.fa.fai --extract-dir split_by_contig
```

## Richer output

```bash
samtools view -h overhangs.bam \
  | teloclip extract \
      --ref-idx ref.fa.fai \
      --extract-dir split_by_contig \
      --include-stats \
      --count-motifs TTAGGG \
      --fuzzy-count \
      --report-stats
```

- `--include-stats` puts MAPQ, clip length and motif counts in each FASTA header
- `--count-motifs` counts the given motifs in the clipped region only
- `--report-stats` prints a summary, including contigs with no overhangs at all

By default the clipped portion of each read is lowercased, so the boundary
between aligned and overhanging sequence is visible at a glance. Pass
`--no-mask` to disable that.

## What to look for

Open a couple of the files and read them.

**Do the clips agree?** Reads at the same contig end should tell the same story.
Several reads whose clipped sequences align to each other is good evidence.
Reads that disagree wildly suggest non-specific alignment — see
[Subtelomeric repeats](../guide/biology.md#subtelomeric-repeats).

**Is the repeat where it should be?** For the canonical vertebrate repeat you
expect `CCCTAA` accumulating at the 5' end and `TTAGGG` at the 3' end. Finding
the same motif at both ends in the same orientation is a warning sign.

**How long are the clips?** A handful of bases is not a telomere. Hundreds of
bases of clean repeat is.

**Is one read carrying the result?** If a single read has a much longer clip
than the rest, check it is not chimeric before letting `extend` choose it.

## Contigs with no overhangs

`--report-stats` lists contigs where nothing was found. That is informative:
a chromosome-level assembly should have overhang support at most chromosome
ends. Ends with nothing may already be complete, or may be where coverage ran
out.

Next: [Extending contigs](extend.md).
