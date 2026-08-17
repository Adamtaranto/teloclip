# CLI reference

Generated from `--help`. Run the commands yourself for the authoritative text.

## teloclip

```text
Usage: teloclip [OPTIONS] [COMMAND] [ARGS]...

  A tool for the recovery of unassembled telomeres from soft-clipped read
  alignments.

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  extend   Extend contigs using overhang analysis from soft-clipped...
  extract  Extract overhanging reads for each end of each reference contig.
  filter   Filter SAM file for clipped alignments containing unassembled...
```

## teloclip filter

Reads SAM on stdin, writes SAM on stdout.

```text
Usage: teloclip filter [OPTIONS] [SAMFILE]

  Filter SAM file for clipped alignments containing unassembled telomeric
  repeats.

Options:
  --ref-idx PATH                  Path to fai index for reference fasta. Index
                                  fasta using `samtools faidx FASTA`
                                  [required]
  --min-clip INTEGER              Require clip to extend past ref contig end
                                  by at least N bases. Default: 1
  --max-break INTEGER             Tolerate max N unaligned bases before contig
                                  end. Default: 50
  --motifs TEXT                   If set keep only reads containing given
                                  motif/s from comma delimited list of
                                  strings. By default also search for reverse
                                  complement of motifs. i.e. TTAGGG,TTAAGGG
                                  will also match CCCTAA,CCCTTAA
  --no-rev                        If set do NOT search for reverse complement
                                  of specified motifs.
  --keep-secondary                If set, include secondary alignments in
                                  output. Default: Off (exclude secondary
                                  alignments).
  --fuzzy                         If set, tolerate +/- 1 variation in motif
                                  homopolymer runs i.e. TTAGGG ->
                                  T{1,3}AG{2,4}. Default: Off
  -r, --min-repeats INTEGER       Minimum number of sequential pattern matches
                                  required for a hit to be reported. Default:
                                  1
  --min-anchor INTEGER            Minimum number of aligned bases (anchor)
                                  required on the non-clipped portion of the
                                  read. Default: 100
  --match-anywhere                If set, motif match may occur in unclipped
                                  region of reads.
  --log-level [debug|info|warning|error]
                                  Logging level (default: INFO).
  --logfile PATH                  Also write log messages to this file (parent
                                  directories are created).
  --help                          Show this message and exit.
```

## teloclip extract

Reads SAM on stdin, writes per-contig-end FASTA files.

```text
Usage: teloclip extract [OPTIONS] [SAMFILE]

  Extract overhanging reads for each end of each reference contig. Reads are
  always written to output files.

Options:
  --ref-idx PATH                  Path to fai index for reference fasta. Index
                                  fasta using `samtools faidx FASTA`
                                  [required]
  --prefix TEXT                   Use this prefix for output files. Default:
                                  None.
  --extract-dir PATH              Write extracted reads to this directory.
                                  Default: cwd.
  --min-clip INTEGER              Require clip to extend past ref contig end
                                  by at least N bases. Default: 1
  --max-break INTEGER             Tolerate max N unaligned bases before contig
                                  end. Default: 50
  --min-anchor INTEGER            Minimum anchored alignment length required
                                  (default: 100).
  --min-mapq INTEGER              Minimum mapping quality required (default:
                                  0).
  --keep-secondary                If set, include secondary alignments in
                                  output. Default: Off (exclude secondary
                                  alignments).
  --include-stats                 Include mapping quality, clip length, and
                                  motif counts in FASTA headers.
  --count-motifs TEXT             Comma-delimited motif sequences to count in
                                  overhang regions (e.g., "TTAGGG,CCCTAA").
  --fuzzy-count                   Use fuzzy motif matching allowing ±1
                                  character variation when counting motifs.
  --buffer-size INTEGER           Number of sequences to buffer before writing
                                  (default: 1000).
  --output-format [fasta|fastq]   Output format for extracted sequences
                                  (default: fasta).
  --report-stats                  Write extraction statistics to file in
                                  output directory.
  --no-mask-overhangs             Do not convert overhang sequences to
                                  lowercase.
  --log-level [DEBUG|INFO|WARNING|ERROR]
                                  Logging level (default: INFO).
  --logfile PATH                  Also write log messages to this file (parent
                                  directories are created).
  --help                          Show this message and exit.
```

## teloclip extend

Reads an indexed BAM and FASTA from disk; cannot read a stream.

```text
Usage: teloclip extend [OPTIONS] BAM_FILE REFERENCE_FASTA

  Extend contigs using overhang analysis from soft-clipped alignments.

Options:
  --output-fasta PATH             Extended FASTA output file
  --stats-report PATH             Statistics report output file
  --outlier-threshold FLOAT       Modified z-score above which a contig end is
                                  reported as having anomalous overhang
                                  coverage (default: 3.5)
  --min-overhangs INTEGER         Minimum supporting overhangs required
                                  (default: 1)
  --max-homopolymer INTEGER       Maximum homopolymer run length allowed
                                  (default: 500)
  --min-extension INTEGER         Minimum novel bases an overhang must
                                  contribute to be used (default: 1)
  --min-clip INTEGER              Require clip to extend past the contig end
                                  by at least N bases (default: 1)
  --max-break INTEGER             Maximum gap allowed between alignment and
                                  contig end (default: 50)
  --min-anchor INTEGER            Minimum anchor length required for alignment
                                  (default: 100)
  --dry-run                       Report extensions without modifying
                                  sequences
  --count-motifs TEXT             Comma-delimited motif sequences to count in
                                  overhang regions (e.g., "TTAGGG,CCCTAA")
  --fuzzy-count                   Use fuzzy motif matching allowing ±1
                                  character variation when counting motifs
  --prefix TEXT                   Prefix for default output filenames
                                  (default: teloclip_extended)
  --screen-terminal-bases INTEGER
                                  Number of terminal bases to screen for
                                  motifs in original contigs (default: 0,
                                  disabled)
  --exclude-contigs TEXT          Comma-delimited list of contig names to
                                  exclude from extension (e.g.,
                                  "chrM,chrC,scaffold_123")
  --exclude-contigs-file PATH     Text file containing contig names to exclude
                                  (one per line)
  --log-level [debug|info|warning|error]
                                  Logging level (default: INFO).
  --logfile PATH                  Also write log messages to this file (parent
                                  directories are created).
  --html-report PATH              Write a self-contained HTML report showing
                                  every overhang read aligned against the
                                  contig terminus it supports, plus overhang
                                  depth across the assembly.
  --html-max-reads INTEGER        Maximum overhang reads rendered per contig
                                  end in the HTML report (default: 25). Reads
                                  contributing the most sequence are shown
                                  first.
  --overhang-log PATH             Write a TSV describing every accepted
                                  overhang read: contig, end, gap from the
                                  contig terminus, clip length and overhang
                                  length.
  --help                          Show this message and exit.
```
