# 2. Mapping reads

Teloclip reads alignments, so the first step is mapping your raw long reads back
onto the draft assembly.

## Index the assembly

```bash
samtools faidx ref.fa
```

This produces `ref.fa.fai`, which every teloclip sub-command needs to know how
long each contig is.

## Align

=== "PacBio CLR"

    ```bash
    minimap2 -ax map-pb ref.fa reads.fq.gz \
      | samtools sort -o raw.bam
    samtools index raw.bam
    ```

=== "PacBio HiFi"

    ```bash
    minimap2 -ax map-hifi ref.fa reads.fq.gz \
      | samtools sort -o raw.bam
    samtools index raw.bam
    ```

=== "Oxford Nanopore"

    ```bash
    minimap2 -ax map-ont ref.fa reads.fq.gz \
      | samtools sort -o raw.bam
    samtools index raw.bam
    ```

## What matters here

**Use the raw reads.** Teloclip recovers sequence the assembler declined to
place. Corrected or trimmed read sets may have had exactly that sequence removed.

**Do not filter out clipped alignments.** Some pipelines discard heavily clipped
reads as low quality. Those are the reads teloclip needs.

**Soft clips, not hard clips.** Reads whose CIGAR contains `H` are skipped,
because the clipped bases are absent from the record. `minimap2` soft-clips
primary alignments by default. Supplementary alignments are hard-clipped by
convention and are skipped.

**Keep the header.** Teloclip's `filter` and `extract` read SAM text on stdin
and pass the header through. Always use `samtools view -h`.

## Check the mapping

```bash
samtools flagstat raw.bam
```

Before going further, confirm that a reasonable fraction of reads mapped. If
mapping rate is low, teloclip will find little, and the problem is upstream.

You can also check that clipped reads exist at all:

```bash
samtools view raw.bam | awk '$6 ~ /S/' | wc -l
```

Zero here means either your aligner is hard-clipping or your reads do not extend
past the contigs.

Next: [Filtering overhangs](filter.md).
