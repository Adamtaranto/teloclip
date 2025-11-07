[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/teloclip.svg)](https://badge.fury.io/py/teloclip)
[![codecov](https://codecov.io/gh/adamtaranto/teloclip/graph/badge.svg?token=NBS8YPLZDT)](https://codecov.io/gh/adamtaranto/teloclip)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/teloclip/README.html)
[![Downloads](https://pepy.tech/badge/teloclip)](https://pepy.tech/project/teloclip)
[![Docker Image](https://img.shields.io/docker/v/adamtaranto/teloclip?label=docker&color=blue)](https://hub.docker.com/r/adamtaranto/teloclip)
[![Docker Pulls](https://img.shields.io/docker/pulls/adamtaranto/teloclip)](https://hub.docker.com/r/adamtaranto/teloclip)

<br clear="right"/>
<br clear="left"/>

<p align="center">
<img src="https://raw.githubusercontent.com/Adamtaranto/teloclip/main/docs/teloclip_hexlogo.jpg" width="180" height="180" title="teloclip_hex" />
</p>

<h1>Teloclip</h1>
<p>
A tool for the recovery of unassembled telomeres from raw long-reads using soft-clipped read alignments.

</p>

<h3>🎉🧬 New Release v0.3.2: Teloclip now supports automatic telomere extension!! 🧬🎉</h3>

### Table of contents

- [About Teloclip](#about-teloclip)
- [CLI Structure](#cli-structure)
- [Options and Usage](#options-and-usage)
  - [Installation](#installation)
- [Example Usage](#example-usage)
  - [Optional Quality Control](#optional-quality-control)
- [Options](#options)
  - [Main Command](#main-command)
  - [Filter Sub-command Options](#filter-sub-command-options)
  - [Extract Sub-command Options](#extract-sub-command-options)
  - [Extend Sub-command Options](#extend-sub-command-options)
- [Citing Teloclip](#citing-teloclip)
- [Publications using Teloclip](#publications-using-teloclip)
- [Issues](#issues)
- [License](#license)

## About Teloclip

In most eukaryotic species, chromosomes terminate in repetitive [telomeric](https://en.wikipedia.org/wiki/Telomere) sequences. A complete genome assembly should ideally comprise chromosome-level contigs that possess telomeric repeats at each end. However, genome assemblers frequently fail to recover these repetitive features, instead producing contigs that terminate immediately prior to telomeric repeats.

Teloclip is designed to scan raw long-read data for evidence that can be used to restore missing telomeres. It does this by searching alignments of raw long-read data (i.e. Pacbio or ONT reads mapped with Minimap2) for 'clipped' alignments that occur at the ends of draft contigs. A 'clipped' alignment is produced where the _end_ of a read is not part of its best alignment. This can occur when a read extends past the end of an assembled contig.

Information about segments of a read that were aligned or clipped are stored in [SAM formatted](<https://en.wikipedia.org/wiki/SAM_(file_format)>) alignments as a [CIGAR string](https://www.drive5.com/usearch/manual/cigar.html). Teloclip parses these strings to determine if a read has been clipped at one or both ends of a contig.

Optionally, teloclip can screen overhanging reads for telomere-associated motifs (i.e. 'TTAGGG' / 'CCCTAA') and report only those containing a match.

Once candidate telomeric sequences have be detected in alignment overhangs, teloclip can be used to automatically patch the missing sequence onto draft contigs.

Teloclip is based on concepts from Torsten Seemann's excellent tool [samclip](https://github.com/tseemann/samclip). Samclip can be used to remove clipped alignments from a samfile prior to variant calling.

## CLI Structure

Teloclip provides three sub-commands:

- **`teloclip filter`**: Filter SAM/BAM files to identify terminal soft-clipped alignments containing potential telomeric sequences
- **`teloclip extract`**: Extract overhanging reads to separate FASTA files organized by contig and end position
- **`teloclip extend`**: Extend draft contigs using overhang analysis from soft-clipped alignments.

## Options and Usage

### Installation

Teloclip requires Python >= 3.8.

There are 5 options available for installing Teloclip locally:

1. Install from PyPi.
   This or Bioconda will get you the latest stable release.

```bash
pip install teloclip
```

2. Install from Bioconda.

```bash
conda install -c bioconda teloclip
```

3. Pip install directly from this git repository.

This is the best way to ensure you have the latest development version.

```bash
pip install git+https://github.com/Adamtaranto/teloclip.git
```

4. Clone from this repository and install as a local Python package.

Do this if you want to edit the code.

```bash
git clone https://github.com/Adamtaranto/teloclip.git && cd teloclip && pip install -e '.[dev]'
```

5. Use Docker for reproducible containerized environments.

Ideal for pipelines and reproducible workflows. No local Python installation required.

```bash
# Pull the latest image
docker pull adamtaranto/teloclip:latest

# Run teloclip
docker run --rm -v $(pwd):/data adamtaranto/teloclip:latest --version
```

See [DOCKER.md](DOCKER.md) for complete Docker usage guide and [examples/nextflow/](examples/nextflow/) for Nextflow integration.

**Verify installation**

```bash
# Print version number and exit.
teloclip --version
# > teloclip 0.3.2

# Get usage information
teloclip --help
```

## Example Usage

Basic use case:

<p align="center">
<img src="https://raw.githubusercontent.com/Adamtaranto/teloclip/main/docs/teloclip_example_graphic.png" title="teloclip_example" />
</p>

First index the reference assembly so teloclip knows where each contig ends.

```bash
# Create index of reference fasta
samtools faidx ref.fa
```

Next align your raw long reads to the reference fasta.

```bash
minimap2 -ax map-pb ref.fa pacbio_reads.fq.gz > in.sam
```

**Loading alignments from file**

Next you will need to provide alignment records to teloclip in SAM format. These can be read directly from a SAM file like this:

```bash
# Option 1: Read alignment input from sam file and write overhang-reads to stdout
teloclip filter --ref-idx ref.fa.fai in.sam

# Option 2: Read alignment input from stdin and write stdout to file
teloclip filter --ref-idx ref.fa.fai < in.sam > overhangs.sam
```

Alternatively, you can read and write alignment records from BAM files.

BAM files are binary SAM files, they contain all the same information but take up much less storage space.

You can use BAM files with teloclip like this:

```bash
# Read alignments from bam file, pipe sam lines to teloclip, sort overhang-read alignments and write to bam file
samtools view -h in.bam | teloclip filter --ref-idx ref.fa.fai | samtools sort > overhangs.bam
```

**Streaming alignments from Minimap**

You can also stream SAM records directly from the aligner to save disk space.

```bash
# Map PacBio long-reads to ref assembly,
# return alignments clipped at contig ends,
# write to sorted bam.
minimap2 -ax map-pb ref.fa pacbio_reads.fq.gz | teloclip filter --ref-idx ref.fa.fai | samtools sort > overhangs.bam
```

**Report clipped alignments containing target motifs**

`teloclip filter` has the option to report only overhanging reads that contain a known telomeric repeat sequence.

```bash
# Report alignments which are clipped at a contig end
# AND contain >=1 copy of the telomeric repeat "TTAGGG" (or its reverse complement "CCCTAA") in the clipped region.
samtools view -h in.bam | teloclip filter --ref-idx ref.fa.fai --motifs TTAGGG | samtools sort > overhangs.bam

# To change the minimum number of consecutive motif repeats required for a match, set "--min-repeats". This example will require one instance of "TTAGGGTTAGGGTTAGGG" in the overhang.
samtools view -h in.bam | teloclip filter --ref-idx ref.fa.fai --motifs TTAGGG --min-repeats 3 | samtools sort > out.bam
```

**Matching noisy target motifs**

Raw long-reads can contain errors in the length of homopolymer tracks. If the `--fuzzy` option is set, motifs will be converted to regex patterns that allow the number of repeated bases to vary by +/- 1.
i.e. "TTAGGG" -> "T{1,3}AG{2,4}". This pattern will match TTAGG TTAGGGG TAGG TTTAGGG etc.

To reduce off target matching you can increase the minimum required number of sequential motif matches with "--min-repeats".

```bash
samtools view -h in.bam | teloclip filter --ref-idx ref.fa.fai --fuzzy --motifs TTAGGG --min-repeats 4 | samtools sort > overhangs.bam
```

**Extract clipped reads**

`teloclip extract` will write overhanging reads to separate fasta files for each reference contig end. The clipped region of each read is masked as lowercase in output fasta files.

You can inspect these reads and select candidates to manually extend contig ends.

```bash
# Find soft-clipped alignments containing motif 'TTAGGG' that overhang contig ends, write to sorted bam.
samtools view -h in.bam | teloclip filter --ref-idx ref.fa.fai --motifs TTAGGG | samtools sort > sorted_overhangs.bam

# Extract overhang reads and write to separate fasta files for each reference contig end.
# Adds overhang stats to fasta header and writes overhang region in lowercase.
# Note: Use sorted input to make processing more efficient.
samtools view -h sorted_overhangs.bam | teloclip extract --ref-idx ref.fa.fai --extract-dir split_overhangs_by_contig --include-stats --count-motifs TTAGGG --report-stats
```

**Automatically extend missing telomeres**

Use the `teloclip extend` tool to automatically extend contigs with missing telomeic sequences from overhang-reads identified with `teloclip filter`.

Before using overhangs identified by Teloclip to extend contigs you should inspect the alignments in a genome browser that displays information about clipped reads, such as [IGV](https://github.com/igvteam/igv).

Check for conflicting soft-clipped sequences. These indicate non-specific read alignments. You may need to tighten your alignment criteria or manually remove low-confidence alignments.

Note: Circular genomes (i.e. mitochondria, chloroplasts, and nitroplasts) will always yield soft-clipped overhangs and should not be extended. We attempt to exclude these with `--exclude-outliers` which will skip contigs with unusually high overhang depths. You can explicitly exclude known circular contigs by providing names to `--exclude-contigs`.

```bash
# Create required indices (one-time setup)
samtools faidx ref.fa

# Convert SAM -> BAM, sort, and write sorted BAM
samtools view -bS overhangs.sam | samtools sort -o overhangs.sorted.bam

# Index the sorted BAM for fast access
samtools index overhangs.sorted.bam

# Use `--dry-run` option to report proposed changes without applying them.
teloclip extend overhangs.sorted.bam ref.fa \
  --output-fasta extended.fasta \
  --stats-report extension_report.txt \
  --count-motifs TTAGGG \
  --screen-terminal-bases 1000 \
  --exclude-contigs ctg_007_mitochondrial
  --dry-run
```

After manually extending contigs the revised assembly should be re-polished using available long and short read data to correct indels present in the raw long-reads.

The final telomere-extended assembly should be re-polished using available long and short read data to correct indels (i.e. with `Medaka` and `Pypolca`) in the raw long-read extensions.

### Optional Quality Control

**Additional filters**

Users may wish to exclude reads below a minimum mapping quality score to reduce the risk of incorrect alignments.

Similarly, multi-mapping reads will generate secondary alignments. To exclude non-specific aligments you can pre-filtering with `samtools view`. You can [decode sam flags here](https://broadinstitute.github.io/picard/explain-flags.html).

Note: As of version teloclip v0.3.0, `filter` and `extract` will exclude secondary alignments by default.

```bash
# Use samtools to filter reads below a MAPQ 30
samtools view -h -q 30 input.sam | teloclip filter --ref-idx ref.fa.fai > min_mapq_30.sam

# Exclude secondary alignments by filtering with samtools
# Note: Secondary alignments are filtered by default in teloclip >=v0.3.0, use '--keep-secondary' to keep.
samtools view -h -F 0x100 input.sam | teloclip filter --ref-idx ref.fa.fai > no_secondary.sam
```

## Options

The main `teloclip` command provides global options and sub-commands for specific operations.

### Main Command

Run `teloclip --help` to view the main command options:

```code
Usage: teloclip [OPTIONS] COMMAND [ARGS]...

  A tool for the recovery of unassembled telomeres from soft-clipped read
  alignments.

Options:
  -v, --verbose                   Enable verbose logging
  -q, --quiet                     Suppress all but error messages
  --log-level [DEBUG|INFO|WARNING|ERROR]
                                  Set specific log level
  --version                       Show the version and exit.
  --help                          Show this message and exit.

Commands:
  extend   Extend contigs using overhang analysis from soft-clipped...
  extract  Extract overhanging reads for each end of each reference contig.
  filter   Filter SAM file for clipped alignments containing unassembled...
```

### Filter Sub-command Options

Run `teloclip filter --help` to view the filter command options:

```code
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
  --log-level [DEBUG|INFO|WARNING|ERROR]
                                  Logging level (default: INFO).
  --help                          Show this message and exit.
```

### Extract Sub-command Options

Run `teloclip extract --help` to view the extract command options:

```code
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
  --help                          Show this message and exit.
```

### Extend sub-command options

Run `teloclip extend --help` to view the extract command options:

```code
Usage: teloclip extend [OPTIONS] BAM_FILE REFERENCE_FASTA

  Extend contigs using overhang analysis from soft-clipped alignments.

Options:
  --output-fasta PATH             Extended FASTA output file
  --stats-report PATH             Statistics report output file
  --exclude-outliers              Exclude outlier contigs from extension
  --outlier-threshold FLOAT       Z-score threshold for outlier detection
                                  (default: 2.0)
  --min-overhangs INTEGER         Minimum supporting overhangs required
                                  (default: 1)
  --max-homopolymer INTEGER       Maximum homopolymer run length allowed
                                  (default: 500)
  --min-extension INTEGER         Minimum overhang length for extension
                                  (default: 1)
  --max-break INTEGER             Maximum gap allowed between alignment and
                                  contig end (default: 10)
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
  --log-level [DEBUG|INFO|WARNING|ERROR]
                                  Logging level (default: INFO).
  --help                          Show this message and exit.
```

## Citing Teloclip

If you use Teloclip in your work please cite this git repo directly and note the release version you used.

## Publications using Teloclip

Teloclip has been used to recover and extend telomeric sequences in a wide variety of taxa, including Algae, Plants, Insects, and Fungi.

- Deng, Y., Zhou, P., Li, F., Wang, J., Xie, K., Liang, H., Wang, C., Liu, B., Zhu, Z., Zhou, W. and Dun, B., **2024**. A complete assembly of the sorghum BTx623 reference genome. Plant Communications, 5(6). 🌾

- He, W., Hu, D., Guo, M., Nie, B., Zhang, G., Jia, Y., Hou, Z., Shu, S., Shao, Y., Simonsen, H.T. and Twamley, A., **2025**. The telomere‐to‐telomere genome of Sanicula chinensis unveils genetic underpinnings of low furanocoumarin diversity and content in one basal lineage of Apiaceae. The Plant Journal, 123(1), p.e70311. 🌱

- Jaiswal, R.K., Garibo Domingo, T., Grunchec, H., Singh, K., Pirooznia, M., Elhaik, E. and Cohn, M., **2025**. Subtelomeric elements provide stability to short telomeres in telomerase-negative cells of the budding yeast Naumovozyma castellii. Current Genetics, 71(1), p.19. 🍄

- Liu, Y., Chen, Y., Ren, Z. et al. Two haplotype-resolved telomere-to-telomere genome assemblies of Xanthoceras sorbifolium. Sci Data 12, 791 (**2025**). 🌿

- Loos, A., Doykova, E., Qian, J., Kümmel, F., Ibrahim, H., Kiss, L., Panstruga, R. and Kusch, S., **2025**. Saprotrophic Arachnopeziza Species as New Resources to Study the Obligate Biotrophic Lifestyle of Powdery Mildew Fungi. Molecular Ecology Resources, p.e70045. 🍄

- Oberti, H., Sessa, L., Oliveira‐Rizzo, C., Di Paolo, A., Sanchez‐Vallet, A., Seidl, M.F. and Abreo, E., 2025. Novel genomic features in entomopathogenic fungus Beauveria bassiana ILB308: accessory genomic regions and putative virulence genes involved in the infection process of soybean pest Piezodorus guildinii. Pest Management Science, 81(4), pp.2323-2336. 🍄

- van Westerhoven, A.C., Mehrabi, R., Talebi, R., Steentjes, M.B., Corcolon, B., Chong, P.A., Kema, G.H. and Seidl, M.F., **2024**. A chromosome-level genome assembly of Zasmidium syzygii isolated from banana leaves. G3: Genes, Genomes, Genetics, 14(3), p.jkad262. 🍄

- Wan, L., Deng, C., Liu, B. et al. Telomere-to-telomere genome assemblies of three silkworm strains with long-term pupal characteristics. Sci Data 12, 501 (**2025**). 🐛

- Wang, Z.F., Yu, E.P., Fu, L., Deng, H.G., Zhu, W.G., Xu, F.X. and Cao, H.L., **2025**. Chromosome-scale assemblies of three Ormosia species: repetitive sequences distribution and structural rearrangement. GigaScience, 14, p.giaf047. 🌿

- Xu, Z., Wang, G., Zhu, X. et al. Genome assembly of two allotetraploid cotton germplasms reveals mechanisms of somatic embryogenesis and enables precise genome editing. Nat Genet 57, 2028–2039 (**2025**). 🌱

- Yang, H.P., Wenzel, M., Hauser, D.A., Nelson, J.M., Xu, X., Eliáš, M. and Li, F.W., **2021**. Monodopsis and Vischeria genomes shed new light on the biology of eustigmatophyte algae. Genome biology and evolution, 13(11), p.evab233. 🦠

## Issues

Submit feedback to the [Issue Tracker](https://github.com/Adamtaranto/teloclip/issues)

## License

Software provided under MIT license.

## Star History

[![Star History
Chart](https://api.star-history.com/svg?repos=adamtaranto/teloclip&type=Date)](https://star-history.com/#adamtaranto/teloclip&Date)
