# 6. Validating the result

An extension is a claim about your genome. This page covers how to check it
before building on it.

## Mechanical checks

These should always pass. If one fails, something is wrong with the run rather
than the biology.

**Every contig is present, in the original order.**

```bash
samtools faidx extended.fasta
diff <(cut -f1 ref.fa.fai) <(cut -f1 extended.fasta.fai) && echo "order preserved"
```

**No contig got shorter.**

```bash
paste <(cut -f1,2 ref.fa.fai) <(cut -f2 extended.fasta.fai) \
  | awk '$3 < $2 { print "SHORTER:", $1, $2, "->", $3; bad=1 }
         END { if (!bad) print "no contig shortened" }'
```

**The report's arithmetic reconciles.** In the Extensions Applied table,
`Original bp + Total +bp` must equal `Final bp` on every row. The `+bp` columns
are net of trimming, so this is a genuine check rather than a tautology.

## Biological checks

**Count the telomeres.** A chromosome-level assembly should approach two
telomeres per chromosome, and no more.

```bash
grep -c "^>" extended.fasta   # contigs
```

Compare **Contig ends extended** in the summary against twice your chromosome
count. Greatly exceeding it means most extended ends are not telomeres.

**Check motif gain is at the ends.** With `--screen-terminal-bases` set, the
report gives per-end motif counts before and after:

| Contig | End | Motif | Window bp | Pre | Post | Gain |
|---|---|---|---:|---:|---:|---:|
| chr1 | left | CCCTAA | 1,162 | 2 | 29 | +27 |
| chr1 | right | TTAGGG | 1,213 | 1 | 35 | +34 |

Gain concentrated at the ends, in the correct orientation per end, is what you
want. `CCCTAA` at the 5' end and `TTAGGG` at the 3' end is the biologically
correct arrangement for the canonical vertebrate repeat; the same motif in the
same orientation at both ends is a warning sign.

**Review flagged contigs.** Anything under **Contigs With Anomalous Overhang
Coverage** needs a decision from you. See
[Biological cases that cause trouble](../guide/biology.md).

**Re-map and inspect.** The strongest check is to map the reads back onto the
extended assembly and look at the new ends in IGV. Extensions supported by many
reads with clean alignments are trustworthy; ones where support drops off
abruptly are not.

```bash
minimap2 -ax map-ont extended.fasta reads.fq.gz | samtools sort -o check.bam
samtools index check.bam
```

## Polish

!!! warning "Do not skip this"

    `extend` splices in sequence from a **single raw read**, which carries that
    read's indel errors. The extension is as error-prone as the read it came
    from.

Re-polish the extended assembly before any downstream analysis:

```bash
# Long-read polishing
# e.g. Medaka (ONT) or NextPolish2

# Short-read polishing, if you have Illumina data
# e.g. Pypolca
```

Polishing after extension also gives an independent check: if reads do not
support the extended sequence, polishing will disagree with it.

## When to stop

Do not extend if:

- The clipped sequences at an end disagree with each other
- Only one read supports the end and it has an unusually long clip
- The contig is circular, or flagged for anomalous coverage and you cannot
  explain why
- The motif is present but the contig end is internal by other evidence

A draft assembly with honest ends is more useful than one with confident wrong
ones.
