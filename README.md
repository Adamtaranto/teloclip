[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI version](https://badge.fury.io/py/teloclip.svg)](https://badge.fury.io/py/teloclip)
[![codecov](https://codecov.io/gh/adamtaranto/teloclip/graph/badge.svg?token=NBS8YPLZDT)](https://codecov.io/gh/adamtaranto/teloclip)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/teloclip.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/teloclip)
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

<h3>🎉🧬 Teloclip now supports automatic telomere extension with <code>teloclip extend</code>!! 🧬🎉</h3>

📖 **[Full documentation](https://adamtaranto.github.io/teloclip/)** — tutorial,
guidance on interpreting results, and CLI reference.

### Table of contents

- [About Teloclip](#about-teloclip)
- [CLI Structure](#cli-structure)
- [Options and Usage](#options-and-usage)
  - [Installation](#installation)
- [Example Usage](#example-usage)
  - [Optional Quality Control](#optional-quality-control)
- [Options](#options)
  - [Filter Sub-command Options](#filter-sub-command-options)
  - [Extend Sub-command Options](#extend-sub-command-options)
- [Citing Teloclip](#citing-teloclip)
- [Publications using Teloclip](#publications-using-teloclip)

## About Teloclip

In most eukaryotic species, chromosomes terminate in repetitive [telomeric](https://en.wikipedia.org/wiki/Telomere) sequences. A complete genome assembly should ideally comprise chromosome-level contigs that possess telomeric repeats at each end. However, genome assemblers frequently fail to recover these repetitive features, instead producing contigs that terminate immediately prior to telomeric repeats.

Part of the reason is a sampling artefact: read coverage cannot stay flat all the way to the end of a linear chromosome. Random fragmentation of linear chromosomes followed by filtering of short fragments causes expected coverage to decay smoothly to zero over roughly one read length. The expected coverage at distance `x` from the terminus is `E[min(x + 1, ℓ)] / E[ℓ]` over the size-selected read-length distribution `ℓ`. Assemblers see decaying coverage as lack of evidence and trim contigs short of the telomere.

<img src="https://raw.githubusercontent.com/Adamtaranto/teloclip/main/docs/images/coverage-decay-bootstrap.png" alt="Simulated coverage against distance from a contig end for lognormal 15 plus or minus 3 Kbp fragments with fragments at or below 10 Kbp discarded, at thirty-fold interior depth. Coverage climbs from zero at the terminus to interior depth about 20 Kbp in, with a 95 percent band from 200 simulations around the mean and the analytic expectation tracking it closely." />

Simulated terminal coverage for a PacBio HiFi-like library (lognormal fragment lengths, 15 ± 3 Kbp, fragments ≤10 Kbp discarded) at 30× interior depth. See [Why coverage dies at contig ends](https://adamtaranto.github.io/teloclip/guide/coverage-decay/) for the full model.

Teloclip is designed to scan raw long-read data for evidence that can be used to restore missing telomeres. It does this by searching alignments of raw long-read data (i.e. Pacbio or ONT reads mapped with Minimap2) for 'clipped' alignments that occur at the ends of draft contigs. A 'clipped' alignment is produced where the _end_ of a read is not part of its best alignment. This can occur when a read extends past the end of an assembled contig.

Information about segments of a read that were aligned or clipped are stored in [SAM formatted](<https://en.wikipedia.org/wiki/SAM_(file_format)>) alignments as a [CIGAR string](https://www.drive5.com/usearch/manual/cigar.html). Teloclip parses these strings to determine if a read has been clipped at one or both ends of a contig.

Optionally, teloclip can screen overhanging reads for telomere-associated motifs (i.e. 'TTAGGG' / 'CCCTAA') and report only those containing a match.

Once candidate telomeric sequences have be detected in alignment overhangs, teloclip can be used to automatically patch the missing sequence onto draft contigs.

> Teloclip is based on concepts from Torsten Seemann's excellent tool [samclip](https://github.com/tseemann/samclip). Samclip can be used to remove clipped alignments from a samfile prior to variant calling.

## CLI Structure

Teloclip provides three sub-commands:

- **`teloclip filter`**: Filter SAM/BAM files to identify terminal soft-clipped alignments containing potential telomeric sequences
- **`teloclip extract`**: Extract overhanging reads to separate FASTA files organized by contig and end position
- **`teloclip extend`**: Extend draft contigs using overhang analysis from soft-clipped alignments.

## Options and Usage

### Installation

Teloclip requires Python >= 3.8.

#### Local Install

Teloclip is typically used alongside `minimap2` and `samtools`.
Use the `environment.yml` file in this repo to create a conda environment with companion tools and dependencies.

```bash
conda env create -f environment.yml
conda activate teloclip
```
Install `Teloclip` into the conda env with:

A. Install from Bioconda.

```bash
conda install -c bioconda teloclip
```

B. Pip install directly from this git repository.

This is the best way to ensure you have the latest development version.

```bash
pip install git+https://github.com/Adamtaranto/teloclip.git
```

**Verify installation**

```bash
# Print version number and exit.
teloclip --version
# > teloclip, version 0.4.1

# Get usage information
teloclip --help
```

#### Run Docker Container

Use Docker for reproducible containerized environments.

Ideal for pipelines and reproducible workflows. No local Python installation required.

```bash
# Pull the latest image
docker pull adamtaranto/teloclip:latest

# Run teloclip
docker run --rm -v $(pwd):/data adamtaranto/teloclip:latest --version
```

See [DOCKER.md](DOCKER.md) for complete Docker usage guide and [examples/nextflow/](examples/nextflow/) for Nextflow integration.

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
# Read alignment input from stdin and write stdout to file
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

**Automatically extend missing telomeres**

Use the `teloclip extend` tool to automatically extend contigs with missing telomeic sequences from overhang-reads identified with `teloclip filter`.

Before using overhangs identified by Teloclip to extend contigs you should inspect the alignments in a genome browser that displays information about clipped reads, such as [IGV](https://github.com/igvteam/igv).

Check for conflicting soft-clipped sequences. These indicate non-specific read alignments. You may need to tighten your alignment criteria or manually remove low-confidence alignments.

Note: Circular genomes (i.e. mitochondria, chloroplasts, and nitroplasts) will always yield soft-clipped overhangs and should not be extended. Exclude them by name with `--exclude-contigs` (or list them one per line and pass `--exclude-contigs-file`).

Teloclip also reports contigs whose overhang coverage is far above the rest of the assembly, which is the usual signature of a collapsed repeat, an rDNA array, or an organellar contig attracting reads from elsewhere. These appear in the stats report under **Contigs With Anomalous Overhang Coverage**. They are *not* excluded automatically: whether extension is appropriate is a judgement about your assembly. Review them and re-run with `--exclude-contigs` if you agree.

Teloclip flags overhang *length* separately from depth. An end that is anomalous on both is the signature of a collapsed array at the terminus; anomalous depth alone more often means an organellar contig or a repeat pulling in reads from elsewhere, and anomalous length alone is often a genuine long telomere that the assembly stopped short of — exactly the case extension exists to serve.

> `--exclude-outliers` has been removed. It previously dropped such contigs silently and the exclusions were never reported anywhere; it was then accepted-but-ignored for a release. Use `--exclude-contigs`.

```bash
# Create required indices (one-time setup)
samtools faidx ref.fa

# Convert SAM -> BAM, sort, and write sorted BAM
samtools view -bS overhangs.sam | samtools sort -o overhangs.sorted.bam

# Index the sorted BAM for fast access
samtools index overhangs.sorted.bam

# Use `--dry-run` to report proposed changes without applying them.
teloclip extend overhangs.sorted.bam ref.fa \
  --output-fasta extended.fasta \
  --stats-report extension_report.md \
  --count-motifs TTAGGG \
  --screen-terminal-bases 1000 \
  --exclude-contigs ctg_007_mitochondrial \
  --dry-run

# Record every overhang read that was considered, and keep a run log.
teloclip extend overhangs.sorted.bam ref.fa \
  --output-fasta extended.fasta \
  --stats-report extension_report.md \
  --overhang-log overhangs.tsv \
  --logfile teloclip_extend.log
```

`--stats-report` accepts a path, or `-` to write the report to stdout. If
omitted it goes to stderr, interleaved with the log.

Extension amounts in the report are **net**: the overhang grafted onto an end,
less any contig bases trimmed to make room for it. So
`Original bp + Total +bp = Final bp` for every row, and you can check the report
against the sequences it describes.

After manually extending contigs the revised assembly should be re-polished using available long and short read data to correct indels present in the raw long-reads.

The final telomere-extended assembly should be re-polished using available long and short read data to correct indels (i.e. with `NextPolish2` and `Pypolca`) in the raw long-read extensions.

### Reading the HTML report

`--html-report` writes a single self-contained file — no CDN, no web fonts, no
separate stylesheet — so it survives being emailed, copied onto a cluster, or
opened years later. The figures below come from the demo assembly built by
`scripts/testdata/generate_demo_assembly.py`, so you can reproduce them
yourself. See [Reading the report](https://adamtaranto.github.io/teloclip/guide/reports/)
for the full walkthrough.

**How much evidence supports each contig end?**

<img src="https://raw.githubusercontent.com/Adamtaranto/teloclip/main/docs/images/report-overhang-depth.png" alt="Strip plot of overhang read depth for each contig end, with three anomalous ends ringed and labelled" />

One point per contig end, blue for left and orange for right, against a median
reference line. A healthy assembly is a flat band near the median. Here
`rdna_plasmid` and `chr7_rdna_array` sit far above it and are ringed as
anomalous.

**How far past the end do those reads reach?**

<img src="https://raw.githubusercontent.com/Adamtaranto/teloclip/main/docs/images/report-overhang-length.png" alt="Split violin plot of overhang length distributions per contig" />

Depth alone will not tell you how much sequence an extension can recover. Each
contig is one shape split at its centre line — left end on the left half, right
end on the right — with the median tick and interquartile range drawn over the
distribution. Eight contigs carry ordinary few-hundred-base overhangs; two run
into kilobases.

**Which kind of end is it?**

<img src="https://raw.githubusercontent.com/Adamtaranto/teloclip/main/docs/images/report-depth-vs-length.png" alt="Scatter plot of overhang depth against median overhang length, with each corner labelled" />

Plotting the two measures against each other separates cases that either one
alone would conflate. Each corner means something different:

| Position | Reading |
|---|---|
| Top right — deep **and** long | A collapsed repeat or rDNA array at the terminus (`chr7_rdna_array`). Extending from a single read is rarely meaningful. |
| Top left — deep but short | Reads drawn in from elsewhere, as a high-copy circular element does (`rdna_plasmid`). |
| Bottom right — long but shallow | A telomere the assembly stopped short of (`chr8_long_telomere`). This is the case extension exists to serve. |
| Bottom left — shallow and short | Little evidence either way. Not an anomaly, just a quiet end. |

Clicking any mark selects that contig in all three charts at once, so an end
that looks unremarkable in one view can be followed to where it is not.

**Should I believe this extension?**

<img src="https://raw.githubusercontent.com/Adamtaranto/teloclip/main/docs/images/report-alignments.png" alt="Per-read alignment panel showing overhang reads laid out against a contig terminus" />

Every supporting read laid out against the contig terminus, marked by the
vertical rule. Grey is the anchored portion of the read, colour is the soft
clip, and highlighted bases are motif matches. This is what lets you judge an
extension rather than take it on trust.


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
  --log-level [debug|info|warning|error]
                                  Logging level (default: INFO).
  --logfile PATH                  Also write log messages to this file (parent
                                  directories are created).
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

## Citing Teloclip

If you use Teloclip in your work please cite this git repo directly and note the release version you used.

## Publications using Teloclip

Teloclip has been used to recover and extend telomeric sequences in a wide variety of taxa, including Algae, Plants, Insects, and Fungi.

- Li, J., Chen, Z., Li, K., Tan, J., Sun, J., Deng, X.W., Park, Y., He, H., Deng, Y. and Zhang, X., **2026**. Telomere-to-telomere genome assembly and a mutant library empower functional genomics and genetic improvement in _Cucurbita moschata_. Plant Communications, 7(5). 🌱

- Liu, Y., Zhao, L., Zhang, J., Ju, Q., Fan, X., Li, Z., Zhang, X., Liang, X., Ge, F. and Chen, J., **2026**. Gapless genome assembly and evolutionary analysis of Cnidium monnieri (Apiaceae). Genomics Communications, 3(1). 🌱

- He, W., Hu, D., Guo, M., Nie, B., Zhang, G., Jia, Y., Hou, Z., Shu, S., Shao, Y., Simonsen, H.T. and Twamley, A., **2025**. The telomere‐to‐telomere genome of _Sanicula chinensis_ unveils genetic underpinnings of low furanocoumarin diversity and content in one basal lineage of Apiaceae. The Plant Journal, 123(1), p.e70311. 🌱

- Jaiswal, R.K., Garibo Domingo, T., Grunchec, H., Singh, K., Pirooznia, M., Elhaik, E. and Cohn, M., **2025**. Subtelomeric elements provide stability to short telomeres in telomerase-negative cells of the budding yeast _Naumovozyma castellii_. Current Genetics, 71(1), p.19. 🍄

- Liu, Y., Chen, Y., Ren, Z. et al. Two haplotype-resolved telomere-to-telomere genome assemblies of _Xanthoceras sorbifolium_. Sci Data 12, 791 (**2025**). 🌿

- Loos, A., Doykova, E., Qian, J., Kümmel, F., Ibrahim, H., Kiss, L., Panstruga, R. and Kusch, S., **2025**. Saprotrophic _Arachnopeziza_ Species as New Resources to Study the Obligate Biotrophic Lifestyle of Powdery Mildew Fungi. Molecular Ecology Resources, p.e70045. 🍄

- Oberti, H., Sessa, L., Oliveira‐Rizzo, C., Di Paolo, A., Sanchez‐Vallet, A., Seidl, M.F. and Abreo, E., **2025**. Novel genomic features in entomopathogenic fungus _Beauveria bassiana_ ILB308: accessory genomic regions and putative virulence genes involved in the infection process of soybean pest _Piezodorus guildinii_. Pest Management Science, 81(4), pp.2323-2336. 🍄

- Wan, L., Deng, C., Liu, B. et al. Telomere-to-telomere genome assemblies of three silkworm strains with long-term pupal characteristics. Sci Data 12, 501 (**2025**). 🐛

- Wang, Z.F., Yu, E.P., Fu, L., Deng, H.G., Zhu, W.G., Xu, F.X. and Cao, H.L., **2025**. Chromosome-scale assemblies of three _Ormosia_ species: repetitive sequences distribution and structural rearrangement. GigaScience, 14, p.giaf047. 🌿

- Xu, Z., Wang, G., Zhu, X. et al. Genome assembly of two allotetraploid cotton germplasms reveals mechanisms of somatic embryogenesis and enables precise genome editing. Nat Genet 57, 2028–2039 (**2025**). 🌱

- Deng, Y., Zhou, P., Li, F., Wang, J., Xie, K., Liang, H., Wang, C., Liu, B., Zhu, Z., Zhou, W. and Dun, B., **2024**. A complete assembly of the sorghum BTx623 reference genome. Plant Communications, 5(6). 🌾

- van Westerhoven, A.C., Mehrabi, R., Talebi, R., Steentjes, M.B., Corcolon, B., Chong, P.A., Kema, G.H. and Seidl, M.F., **2024**. A chromosome-level genome assembly of _Zasmidium syzygii_ isolated from banana leaves. G3: Genes, Genomes, Genetics, 14(3), p.jkad262. 🍄

- Yang, H.P., Wenzel, M., Hauser, D.A., Nelson, J.M., Xu, X., Eliáš, M. and Li, F.W., **2021**. _Monodopsis_ and _Vischeria_ genomes shed new light on the biology of eustigmatophyte algae. Genome biology and evolution, 13(11), p.evab233. 🦠

## Star History

<a href="https://www.star-history.com/?repos=Adamtaranto%2Fteloclip&type=date&legend=top-left">
 <picture>
   <source media="(prefers-color-scheme: dark)" srcset="https://api.star-history.com/chart?repos=Adamtaranto/teloclip&type=date&theme=dark&legend=top-left&sealed_token=jCjxdn3HLLirH-ScvCrYTN-0ifrY4W7QfxFw9K0PsgyDTl_ibIlzfv5kFzeGd_7mzLL14Ftmm6SeBzyKBXHzdMOsjkyJERa1irqLV-QpsqZ6FlcKyytaZOLLbM72VECXt78_7k4-SXz6eh3F-nRAyMZMP6mBOXZJ4TyRHAAae4-w4s8SZ1Bx_F7-kgN2" />
   <source media="(prefers-color-scheme: light)" srcset="https://api.star-history.com/chart?repos=Adamtaranto/teloclip&type=date&legend=top-left&sealed_token=jCjxdn3HLLirH-ScvCrYTN-0ifrY4W7QfxFw9K0PsgyDTl_ibIlzfv5kFzeGd_7mzLL14Ftmm6SeBzyKBXHzdMOsjkyJERa1irqLV-QpsqZ6FlcKyytaZOLLbM72VECXt78_7k4-SXz6eh3F-nRAyMZMP6mBOXZJ4TyRHAAae4-w4s8SZ1Bx_F7-kgN2" />
   <img alt="Star History Chart" src="https://api.star-history.com/chart?repos=Adamtaranto/teloclip&type=date&legend=top-left&sealed_token=jCjxdn3HLLirH-ScvCrYTN-0ifrY4W7QfxFw9K0PsgyDTl_ibIlzfv5kFzeGd_7mzLL14Ftmm6SeBzyKBXHzdMOsjkyJERa1irqLV-QpsqZ6FlcKyytaZOLLbM72VECXt78_7k4-SXz6eh3F-nRAyMZMP6mBOXZJ4TyRHAAae4-w4s8SZ1Bx_F7-kgN2" />
 </picture>
</a>
