# 5. Extending contigs

`teloclip extend` selects the best overhang at each contig end and splices it on,
writing a new assembly.

Unlike `filter` and `extract`, it needs an **indexed BAM** and an **indexed
FASTA** on disk. It cannot read a stream, so `filter | extend` does not work.

## Prepare the inputs

```bash
samtools sort overhangs.sam -o overhangs.bam   # if you have SAM
samtools index overhangs.bam
samtools faidx ref.fa
```

## Dry run first

```bash
teloclip extend overhangs.bam ref.fa \
  --dry-run \
  --stats-report - \
  --count-motifs TTAGGG \
  --screen-terminal-bases 1000
```

`--dry-run` reports exactly what would change without writing sequence.
`--stats-report -` sends the report to stdout; give it a path to write a file, or
omit it and the report goes to stderr interleaved with the log.

Dry-run figures account for trimming, so the predicted final lengths match what
a real run produces.

## Apply it

```bash
teloclip extend overhangs.bam ref.fa \
  --output-fasta extended.fasta \
  --stats-report extension_report.md \
  --count-motifs TTAGGG \
  --screen-terminal-bases 1000 \
  --max-homopolymer 100 \
  --exclude-contigs chrM,chrC \
  --overhang-log overhangs.tsv \
  --logfile extend.log
```

Output contig order matches the input assembly, including contigs that were
excluded or not extended.

## How the overhang is chosen

At each end, candidates are ranked by **net gain**: the clipped bases lying past
the contig terminus, after accounting for any contig bases that must be trimmed
to graft the clip on.

Ranking on net gain rather than raw clip length matters. A 200 bp clip beginning
190 bases inside the contig contributes 10 novel bases; a 150 bp clip flush with
the terminus contributes all 150.

Where two candidates offer comparable gain, the one that trims **less** wins.
Trimming discards polished assembly consensus and replaces it with one raw
read's version of that region, so gaining a single extra base does not justify
discarding thirteen bases of consensus.

Candidates containing a homopolymer run at or above `--max-homopolymer` are
rejected, and the next best is tried.

## Key options

| Option | Default | What it controls |
|---|---|---|
| `--min-extension` | 1 | Novel bases an overhang must contribute to be used |
| `--min-overhangs` | 1 | Supporting reads required before a contig is considered |
| `--max-homopolymer` | 500 | Reject clips containing a run this long |
| `--screen-terminal-bases` | 0 | Terminal window scanned for motifs before and after |
| `--exclude-contigs` | — | Comma-separated names to leave untouched |
| `--overhang-log` | — | TSV of every accepted overhang |

`--min-overhangs 3` is a cheap way to require corroboration.

Unknown `--exclude-contigs` names are an error, not a warning — a typo will not
silently leave a contig in.

!!! note "`--exclude-outliers` is deprecated"

    It is accepted but does nothing. Contigs with anomalous overhang coverage
    are now reported for review instead of being dropped silently. See
    [Reading the report](../guide/reports.md).

## The per-overhang log

`--overhang-log` writes one row per accepted overhang:

```
contig  contig_length  end  read   aln_start  aln_end  gap_from_end  clip_length  overhang_length  anchor_length
chr1    24518201       L    m54..  3          18402    2             164          162              18400
```

This is the raw material behind every decision: which reads were available, how
far each fell short of the terminus, and how much novel sequence each offered.
Useful when a result surprises you.

Next: [Validating the result](validate.md).
