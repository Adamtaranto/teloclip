# teloclip

Find soft-clipped alignments containing unassembled telomeric repeats.

# Table of contents

* [About teloclip](#about-teloclip)
* [Options and usage](#options-and-usage)
    * [Installation](#installation)
    * [Usage](usage)
        * [Examples](examples)
    * [Options](teloclip-options)
* [Issues](issues)
* [License](#license)


## About teloclip

In most eukaryotic species, chromosomes terminate in repetitive [telomeric](https://en.wikipedia.org/wiki/Telomere) 
sequences. A complete genome assembly should ideally comprise chromosome-level contigs that possess telomeric 
repeats at each end. However, genome assemblers frequently fail to recover these repetitive features, instead 
producing contigs that terminate immediately prior to their location.

Teloclip is designed to recover long-reads that can be used to extend draft contigs and resolve missing telomeres 
(short-read alignments may also be processed with teloclip). It does this by searching alignments of raw 
long-read data (i.e. Pacbio or ONP reads mapped with Minimap2) for 'clipped' alignments that occur at the ends of 
draft contigs. A 'clipped' alignment is produced where the *end* of a read is not part of its best alignment. 
This can occur when a read extends past the end of an assembled contig.

Information about segments of a read that were aligned or clipped are stored in [SAM formatted](https://en.wikipedia.org/wiki/SAM_(file_format))
alignments as a [CIGAR string](https://www.drive5.com/usearch/manual/cigar.html). Teloclip parses these strings 
to determine if a read has been clipped at one or both ends of a contig. 

Optionally, teloclip can screen overhanging reads for telomere-associated motifs (i.e. 'TTAGGG' / 'CCCTAA')
and report only those containing a match.

Teloclip is based on concepts from Torsten Seemann's excellent tool [samclip](https://github.com/tseemann/samclip).
Samclip can be used to remove clipped alignments from a samfile prior to variant calling.


## Installation

Clone from this repository and install as a local Python package.

```bash
% git clone https://github.com/Adamtaranto/teloclip.git && cd teloclip && pip install -e .
```

Test installation.

```bash
# Print version number and exit.
% teloclip --version
teloclip 0.0.1

# Get usage information
% teloclip --help
```

## Usage

![teloclip_example](docs/teloclip_example_graphic.png)

*Additional filters*  

  - Consider pre-filtering alignments with "samtools view" to remove non-primary 
    / low quality alignments.
  - Users may wish to exclude reads below a minimum length or read quality score 
    to reduce the risk of incorrect alignments.

*Extending contigs*  

  - Before using terminal alignments identified by teloclip to extend contigs, 
    inspect alignments in a genome browser that displays information about clipped 
    reads, such as [IGV](https://github.com/igvteam/igv). Check for conflicting 
    clipped sequences.
  - After manually extending contigs the revised assembly should be re-polished 
    using available long and short read data to correct indels present in the raw 
    long-reads.
  - Validate the final assembly by re-mapping long-read data and checking for 
    alignments that extend into revised contig ends.

### Examples

teloclip requires an indexed reference fasta
```
# Create index of reference fasta
% samtools faidx ref.fa
```

Read alignments from SAM file
```
# Read input from file and write output to stout
% teloclip --ref ref.fa.fai in.sam

# Read input from stdin and write stdout to file
% teloclip --ref ref.fa.fai < in.sam > out.sam

# Filter alignments from BAM file, write sorted output to file
% samtools view -h in.bam | teloclip --ref ref.fa.fai | samtools sort > out.bam
```

Stream SAM records from aligner
```
# Map PacBio long-reads to ref assembly, filter for alignments clipped at contig ends, write to sorted bam
% minimap2 -ax map-pb ref.fa pacbio.fq.gz | teloclip --ref ref.fa.fai | samtools sort > out.bam 

# Map reads, exclude unmapped reads and non-primary/supplementary alignments. Report clipped reads as sorted bam.
% minimap2 -ax map-pb ref.fa pacbio.fq.gz | samtools view -h -F 0x2308 | teloclip --ref ref.fa.fai | samtools sort > out.bam 

# Map long-reads with MiniMap2 and retain only reads which extend past a cotig end
# AND contain >=1 copy of the telomeric repeat "TTAGGG" (or its reverse complement "CCCTAA")
% minimap2 -ax map-pb ref.fa pacbio.fq.gz | teloclip --ref ref.fa.fai --motifs TTAGGG | samtools sort > out.bam 

```

## Options

Run `teloclip --help` to view the programs' most commonly used options:

```
Usage: teloclip [-h] --refIdx REFIDX [--minClip MINCLIP] [--maxBreak MAXBREAK]
                [--motifs MOTIFS]
                [samfile]

Required:
 --refIdx REFIDX       Path to fai index for reference fasta. Index fasta using `samtools faidx FASTA`

Positional arguments:
  samfile               Input SAM can be added as the first positional argument after flagged options. 
                          If not set teloclip will read from stdin.

Optional:
  --minClip MINCLIP    Require clip to extend past ref contig end by at least N bases.
                         Default: 1
  --maxBreak MAXBREAK  Tolerate max N unaligned bases at contig ends. 
                         Default: 50
  --motifs MOTIFS      If set keep only reads containing given motif/s from a comma delimited list 
                         of strings. By default also search for reverse complement of motifs. 
                         i.e. TTAGGG,TTAAGGG will also match CCCTAA,CCCTTAA
                         Default: None
  --norev NOREV        If set do NOT search for reverse complement of specified motifs. 
                         Default: False
  --matchAny MATCHANY  If set motif match may occur in unclipped region of alignment.
                         Defaut: False
  --version            Show program's version number and exit.
```

## Issues

Submit feedback to the [Issue Tracker](https://github.com/Adamtaranto/teloclip/issues)

## License

Software provided under MIT license.

## Author

[Adam Taranto](https://github.com/Adamtaranto)
