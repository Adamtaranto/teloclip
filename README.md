[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/teloclip.svg)](https://badge.fury.io/py/teloclip)
[![codecov](https://codecov.io/gh/adamtaranto/teloclip/graph/badge.svg?token=NBS8YPLZDT)](https://codecov.io/gh/adamtaranto/teloclip)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/teloclip/README.html)
[![Downloads](https://pepy.tech/badge/teloclip)](https://pepy.tech/project/teloclip)

<br clear="right"/>
<br clear="left"/>

<p align="center">
<img src="https://raw.githubusercontent.com/Adamtaranto/teloclip/main/docs/teloclip_hexlogo.jpg" width="180" height="180" title="teloclip_hex" />
</p>

<h1>Teloclip</h1>
<p>
A tool for the recovery of unassembled telomeres from raw long-reads using soft-clipped read alignments.
</p>

### Table of contents

- [About Teloclip](#about-teloclip)
- [Options and Usage](#options-and-usage)
  - [Installation](#installation)
  - [Run with Gitpod](#run-with-gitpod)
- [Example Usage](#example-usage)
  - [Optional Quality Control](#optional-quality-control)
  - [Extending contigs](#extending-contigs)
  - [Alternative use cases](#alternative-use-cases)
- [Options](#options)
  - [Teloclip Options](#teloclip-options)
  - [Teloclip-extract Options](#teloclip-extract-options)
- [Citing Teloclip](#citing-teloclip)
- [Publications using Teloclip](#publications-using-teloclip)
- [Issues](#issues)
- [License](#license)

## About Teloclip

In most eukaryotic species, chromosomes terminate in repetitive [telomeric](https://en.wikipedia.org/wiki/Telomere) sequences. A complete genome assembly should ideally comprise chromosome-level contigs that possess telomeric repeats at each end. However, genome assemblers frequently fail to recover these repetitive features, instead producing contigs that terminate immediately prior to their location.

Teloclip is designed to recover long-reads that can be used to extend draft contigs and resolve missing telomeres (short-read alignments may also be processed with teloclip). It does this by searching alignments of raw long-read data (i.e. Pacbio or ONT reads mapped with Minimap2) for 'clipped' alignments that occur at the ends of draft contigs. A 'clipped' alignment is produced where the *end* of a read is not part of its best alignment. This can occur when a read extends past the end of an assembled contig.

Information about segments of a read that were aligned or clipped are stored in [SAM formatted](https://en.wikipedia.org/wiki/SAM_(file_format)) alignments as a [CIGAR string](https://www.drive5.com/usearch/manual/cigar.html). Teloclip parses these strings to determine if a read has been clipped at one or both ends of a contig.

Optionally, teloclip can screen overhanging reads for telomere-associated motifs (i.e. 'TTAGGG' / 'CCCTAA') and report only those containing a match.

Teloclip is based on concepts from Torsten Seemann's excellent tool [samclip](https://github.com/tseemann/samclip). Samclip can be used to remove clipped alignments from a samfile prior to variant calling.

## Options and Usage

### Installation

Teloclip requires Python >= 3.8.

There are 4 options available for installing Teloclip locally:

1) Install from PyPi.
This or Bioconda will get you the latest stable release.

```bash
pip install teloclip
```

2) Install from Bioconda.

```bash
conda install -c bioconda teloclip
```

3) Pip install directly from this git repository.

This is the best way to ensure you have the latest development version.

```bash
pip install git+https://github.com/Adamtaranto/teloclip.git
```

4) Clone from this repository and install as a local Python package.

Do this if you want to edit the code.

```bash
git clone https://github.com/Adamtaranto/teloclip.git && cd teloclip && pip install -e '.[dev]'
```

**Verify installation**

```bash
# Print version number and exit.
teloclip --version
# > teloclip 0.1.1

# Get usage information
teloclip --help
```

## Example Usage

Basic use case:

<p align="center">
<img src="https://raw.githubusercontent.com/Adamtaranto/teloclip/main/docs/teloclip_example_graphic.png" title="teloclip_example" />
</p>

**First index the reference assembly**

```bash
# Create index of reference fasta
samtools faidx ref.fa
```

**Reading alignments from SAM file**

```bash
# Read alignment input from sam file and write overhang-reads to stout
teloclip --ref-idx ref.fa.fai in.sam

# Read alignment input from stdin and write stdout to file
teloclip --ref-idx ref.fa.fai < in.sam > out.sam
```

**Reading and writing BAM alignments**

BAM files are binary sam files, they contain all the same information but take up much less storage space.
You can use bam files with teloclip like this:

```bash
# Read alignments from bam file, pipe sam lines to teloclip, sort overhang-read alignments and wite to bam file
samtools view -h in.bam | teloclip --ref-idx ref.fa.fai | samtools sort > out.bam
```

**Streaming SAM records from aligner**

```bash
# Map PacBio long-reads to ref assembly,
# return alignments clipped at contig ends, 
# write to sorted bam.
minimap2 -ax map-pb ref.fa pacbio_reads.fq.gz | teloclip --ref-idx ref.fa.fai | samtools sort > out.bam 

# Map reads to reference, 
# Exclude non-primary alignments. 
# Return alignments clipped at contig ends,
# write to sorted bam.
minimap2 -ax map-pb ref.fa pacbio_reads.fq.gz | samtools view -h -F 0x100 | teloclip --ref-idx ref.fa.fai | samtools sort > out.bam 
```

**Report clipped alignments containing target motifs**

```bash
# Report alignments which are clipped at a contig end
# AND contain >=1 copy of the telomeric repeat "TTAGGG" (or its reverse complement "CCCTAA") in the clipped region.
samtools view -h in.bam | teloclip --ref-idx ref.fa.fai --motifs TTAGGG | samtools sort > out.bam 

# Report alignments which are clipped at a contig end
# AND contain >=1 copy of the telomeric repeat "TTAGGG" (or its reverse complement "CCCTAA") ANYWHERE in the read.
samtools view -h in.bam | teloclip --ref-idx ref.fa.fai --motifs TTAGGG --match-anywhere | samtools sort > out.bam

# To change the minimum number of consecutive motif repeats required for a match, set "--min-repeats"
samtools view -h in.bam | teloclip --ref-idx ref.fa.fai --motifs TTAGGG --min-repeats 4 | samtools sort > out.bam 

```

**Matching noisy target motifs**

Raw long-reads can contain errors in the length of homopolymer tracks. If the `--fuzzy` option is set, motifs will be converted to regex patterns that allow the number of repeated bases to vary by +/- 1.
i.e. "TTAGGG" -> "T{1,3}AG{2,4}". This pattern will match TTAGG TTAGGGG TAGG TTTAGGG etc.

To reduce off target matching you can increase the minimum required number of sequential motif matches with "--min-repeats".

```bash
samtools view -h in.bam | teloclip --ref-idx ref.fa.fai --fuzzy --motifs TTAGGG --min-repeats 4 | samtools sort > out.bam
```

**Extract clipped reads**

`teloclip-extract` will write overhanging reads to separate fasta files for each reference contig end. The clipped region of each read is masked as lowercase in output fasta files.

Collections of reads that overhang a contig end can be assembled with `miniasm` into a single segment before being used to extend the contig. The final telemere-extended assembly should be polished (i.e. with `Racon` or `Pilon`) to correct errors in the raw long-read extensions.

```bash
# Find clipped alignments containing motif 'TTAGGG' and write reads to separate fasta files for each reference contig end.
samtools view -h in.bam | teloclip --ref-idx ref.fa.fai --motifs TTAGGG | teloclip-extract --ref-idx ref.fa.fai --extract-reads --extract-dir split_overhangs_by_contig
```

### Optional Quality Control

**Additional filters**  

Users may wish to exclude reads below a minimum length or read quality score to reduce the risk of incorrect alignments.

In some cases it may be also be useful to prioritise primary alignments. This can be done by pre-filtering alignments with `samtools view`. You can decode sam flags [here](https://broadinstitute.github.io/picard/explain-flags.html).

```bash
# Exclude secondary alignments.
samtools view -h -F 0x100 in.sam | teloclip --ref-idx ref.fa.fai > noSA.sam 
```

**Pre-corrected Data**  

Some assembly tools, such as [Canu](https://github.com/marbl/canu), preform pre-correction of long-reads through iterative overlapping and correction prior to assembly. Corrected reads are trimmed based on coverage to remove low-confidence ends.

This trimming step can result in loss of distal telomeric sequences and so these reads should **NOT** be used with Teloclip.

However, long-reads that have been error-corrected using Illumina data with tools such as [LoRDEC](https://github.com/lanl001/halc) or [HALC](http://www.atgc-montpellier.fr/lordec/) should be fine.

Generally speaking, raw long-reads will be fine for extending your contigs. Any errors in the extended region can be corrected with a round of polishing with short-read data using [Pilon](https://github.com/broadinstitute/pilon).
  
### Extending contigs

Before using terminal alignments identified by Teloclip to extend contigs you should inspect the alignments in a genome browser that displays information about clipped reads, such as [IGV](https://github.com/igvteam/igv).

Check for conflicting soft-clipped sequences. These indicate non-specific read alignments. You may need to tighten your alignment criteria or manually remove low-confidence alignments.

After manually extending contigs the revised assembly should be re-polished using available long and short read data to correct indels present in the raw long-reads.

Finally, validate the updated assembly by re-mapping long-read data and checking for alignments that extend into revised contig ends.

## Options

### Teloclip Options

Run `teloclip --help` to view the programs' most commonly used options:

```code
usage: teloclip [-h] --ref-idx REF_IDX [--min-clip MIN_CLIP] [--max-break MAX_BREAK]
                [--motifs MOTIFS] [--no-rev] [--fuzzy] [-r MIN_REPEATS]
                [--min-anchor MIN_ANCHOR] [--match-anywhere] [-v]
                [samfile]

Filter SAM file for clipped alignments containing unassembled telomeric repeats.

positional arguments:
  samfile

options:
  -h, --help            show this help message and exit
  --ref-idx REF_IDX     Path to fai index for reference fasta. Index fasta using `samtools
                        faidx FASTA`
  --min-clip MIN_CLIP   Require clip to extend past ref contig end by at least N bases.
  --max-break MAX_BREAK
                        Tolerate max N unaligned bases before contig end.
  --motifs MOTIFS       If set keep only reads containing given motif/s from comma
                        delimited list of strings. By default also search for reverse
                        complement of motifs. i.e. TTAGGG,TTAAGGG will also match
                        CCCTAA,CCCTTAA
  --no-rev              If set do NOT search for reverse complement of specified motifs.
  --fuzzy               If set, tolerate +/- 1 variation in motif homopolymer runs i.e.
                        TTAGGG -> T{1,3}AG{2,4}. Default: Off
  -r MIN_REPEATS, --min-repeats MIN_REPEATS
                        Minimum number of sequential pattern matches required for a hit to
                        be reported. Default: 3
  --min-anchor MIN_ANCHOR
                        Minimum number of aligned bases (anchor) required on the non-
                        clipped portion of the read. Default: 500
  --match-anywhere      If set, motif match may occur in unclipped region of reads.
  -v, --version         show program's version number and exit
```

### Teloclip-extract Options

Run `teloclip-extract --help` to view the programs' most commonly used options:

```code
usage: teloclip-extract [-h] --ref-idx REF_IDX [--prefix PREFIX] [--extract-reads]
                        [--extract-dir EXTRACT_DIR] [--min-clip MIN_CLIP]
                        [--max-break MAX_BREAK] [-v]
                        [samfile]

Extract overhanging reads for each end of each reference contig. Write to fasta.

positional arguments:
  samfile

options:
  -h, --help            show this help message and exit
  --ref-idx REF_IDX     Path to fai index for reference fasta. Index fasta using `samtools
                        faidx FASTA`
  --prefix PREFIX       Use this prefix for output files. Default: None.
  --extract-reads       If set, write overhang reads to fasta by contig.
  --extract-dir EXTRACT_DIR
                        Write extracted reads to this directory. Default: cwd.
  --min-clip MIN_CLIP   Require clip to extend past ref contig end by at least N bases.
  --max-break MAX_BREAK
                        Tolerate max N unaligned bases before contig end.
  -v, --version         show program's version number and exit
```

## Citing Teloclip

If you use Teloclip in your work please cite this git repo directly and note the release version you used.

## Publications using Teloclip

Teloclip has been used to recover and extend telomeric sequences in a wide variety of taxa, including Algae, Plants, Insects, and Fungi

Deng, Y., Zhou, P., Li, F., Wang, J., Xie, K., Liang, H., Wang, C., Liu, B., Zhu, Z., Zhou, W. and Dun, B., 2024. A complete assembly of the sorghum BTx623 reference genome. Plant Communications, 5(6).

He, W., Hu, D., Guo, M., Nie, B., Zhang, G., Jia, Y., Hou, Z., Shu, S., Shao, Y., Simonsen, H.T. and Twamley, A., 2025. The telomere‐to‐telomere genome of Sanicula chinensis unveils genetic underpinnings of low furanocoumarin diversity and content in one basal lineage of Apiaceae. The Plant Journal, 123(1), p.e70311.

Jaiswal, R.K., Garibo Domingo, T., Grunchec, H., Singh, K., Pirooznia, M., Elhaik, E. and Cohn, M., 2025. Subtelomeric elements provide stability to short telomeres in telomerase-negative cells of the budding yeast Naumovozyma castellii. Current Genetics, 71(1), p.19.

Liu, Y., Chen, Y., Ren, Z. et al. Two haplotype-resolved telomere-to-telomere genome assemblies of Xanthoceras sorbifolium. Sci Data 12, 791 (2025).

Loos, A., Doykova, E., Qian, J., Kümmel, F., Ibrahim, H., Kiss, L., Panstruga, R. and Kusch, S., 2025. Saprotrophic Arachnopeziza Species as New Resources to Study the Obligate Biotrophic Lifestyle of Powdery Mildew Fungi. Molecular Ecology Resources, p.e70045.

Loos, A., Doykova, E., Qian, J., Kümmel, F., Ibrahim, H., Kiss, L., Panstruga, R. and Kusch, S., 2025. Resources for molecular studies of unculturable obligate biotrophic fungal plant pathogens using their saprotrophic relatives. bioRxiv, pp.2025-05.

Oberti, H., Sessa, L., Oliveira‐Rizzo, C., Di Paolo, A., Sanchez‐Vallet, A., Seidl, M.F. and Abreo, E., 2025. Novel genomic features in entomopathogenic fungus Beauveria bassiana ILB308: accessory genomic regions and putative virulence genes involved in the infection process of soybean pest Piezodorus guildinii. Pest Management Science, 81(4), pp.2323-2336.

van Westerhoven, A.C., Mehrabi, R., Talebi, R., Steentjes, M.B., Corcolon, B., Chong, P.A., Kema, G.H. and Seidl, M.F., 2024. A chromosome-level genome assembly of Zasmidium syzygii isolated from banana leaves. G3: Genes, Genomes, Genetics, 14(3), p.jkad262.

Wan, L., Deng, C., Liu, B. et al. Telomere-to-telomere genome assemblies of three silkworm strains with long-term pupal characteristics. Sci Data 12, 501 (2025).

Wang, Z.F., Yu, E.P., Fu, L., Deng, H.G., Zhu, W.G., Xu, F.X. and Cao, H.L., 2025. Chromosome-scale assemblies of three Ormosia species: repetitive sequences distribution and structural rearrangement. GigaScience, 14, p.giaf047.

Xu, Z., Wang, G., Zhu, X. et al. Genome assembly of two allotetraploid cotton germplasms reveals mechanisms of somatic embryogenesis and enables precise genome editing. Nat Genet 57, 2028–2039 (2025).

Yang, H.P., Wenzel, M., Hauser, D.A., Nelson, J.M., Xu, X., Eliáš, M. and Li, F.W., 2021. Monodopsis and Vischeria genomes shed new light on the biology of eustigmatophyte algae. Genome biology and evolution, 13(11), p.evab233.

## Issues

Submit feedback to the [Issue Tracker](https://github.com/Adamtaranto/teloclip/issues)

## License

Software provided under MIT license.

## Star History

[![Star History
Chart](https://api.star-history.com/svg?repos=adamtaranto/teloclip&type=Date)](https://star-history.com/#adamtaranto/teloclip&Date)
