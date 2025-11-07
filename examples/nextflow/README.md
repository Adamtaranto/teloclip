# Teloclip Nextflow Examples

This directory contains example Nextflow workflows and nf-core style modules for using teloclip in containerized pipelines.

## Quick Start

### Example Workflow

Run the complete example workflow:

```bash
nextflow run teloclip.nf \
  --bam input.bam \
  --ref reference.fa \
  --outdir results \
  --motifs TTAGGG \
  --min_repeats 1
```

### Using Modules

Import teloclip modules into your own workflow:

```groovy
include { TELOCLIP_FILTER } from './modules/teloclip/filter/main'
include { TELOCLIP_EXTRACT } from './modules/teloclip/extract/main'
include { TELOCLIP_EXTEND } from './modules/teloclip/extend/main'

workflow {
    // Filter for telomeric overhangs
    TELOCLIP_FILTER(sam_ch, fai_ch)

    // Sort and index filtered reads
    SAM_TO_BAM(TELOCLIP_FILTER.out.sam)

    // Optional: Extract sequences for inspection
    TELOCLIP_EXTRACT(SAM_TO_BAM.out, fai_ch)

    // Extend reference genome with telomeric sequences
    TELOCLIP_EXTEND(SAM_TO_BAM.out, ref_ch, fai_ch)
}
```

## Modules

### TELOCLIP_FILTER

Filters SAM/BAM files to identify terminal soft-clipped alignments containing potential telomeric sequences.

**Inputs:**

- `meta`: Sample metadata map
- `sam`: SAM format alignment file
- `fai`: FAI index file for reference genome

**Outputs:**

- `sam`: Filtered SAM file containing overhanging reads
- `versions`: Software versions

**Example:**

```groovy
process {
    withName: TELOCLIP_FILTER {
        ext.args = '--motifs TTAGGG --min-repeats 2 --fuzzy'
    }
}
```

### TELOCLIP_EXTRACT

Extracts overhanging reads to separate FASTA files organized by contig and end position.

**Inputs:**

- `meta`: Sample metadata map
- `bam`: BAM format alignment file
- `fai`: FAI index file for reference genome

**Outputs:**

- `fastas`: FASTA files containing overhang sequences
- `versions`: Software versions

**Example:**

```groovy
process {
    withName: TELOCLIP_EXTRACT {
        ext.args = '--min-length 50'
        ext.prefix = 'sample1'
    }
}
```

### TELOCLIP_EXTEND

Extends draft contigs using overhang analysis from soft-clipped read alignments containing telomeric sequences.

**Inputs:**

- `meta`: Sample metadata map
- `bam`: Sorted and indexed BAM file containing filtered soft-clipped alignments (from TELOCLIP_FILTER)
- `bai`: BAM index file
- `fasta`: Reference genome FASTA file
- `fai`: FASTA index file

**Outputs:**

- `fasta`: Extended reference genome FASTA file with telomeres added
- `report`: Statistics report detailing extension operations
- `versions`: Software versions

**Example:**

```groovy
include { TELOCLIP_EXTEND } from './modules/teloclip/extend/main'

process {
    withName: TELOCLIP_EXTEND {
        ext.args = '--count-motifs TTAGGG --screen-terminal-bases 1000 --exclude-outliers'
    }
}
```

## Container Strategy

Each tool runs in its own container:

- **teloclip**: `adamtaranto/teloclip:latest`
- **samtools**: `quay.io/biocontainers/samtools:1.18--h50ea8bc_1`

This follows Nextflow best practices for:

- Reproducibility
- Version control
- Minimal container sizes
- Flexibility in tool versions

## Complete Pipeline Example

```groovy
workflow TELOMERE_ANALYSIS {
    // 1. Index reference
    INDEX_FASTA(ref_ch)

    // 2. Convert BAM to SAM (samtools container)
    BAM_TO_SAM(bam_ch)

    // 3. Filter for telomeric overhangs (teloclip container)
    TELOCLIP_FILTER(BAM_TO_SAM.out, INDEX_FASTA.out)

    // 4. Convert to sorted and indexed BAM (samtools container)
    SAM_TO_BAM(TELOCLIP_FILTER.out.sam)

    // 5. Extract sequences for inspection (teloclip container, optional)
    TELOCLIP_EXTRACT(SAM_TO_BAM.out, INDEX_FASTA.out)

    // 6. Extend reference genome with telomeric sequences (teloclip container)
    TELOCLIP_EXTEND(SAM_TO_BAM.out, ref_ch, INDEX_FASTA.out)
}
```

## Parameters

### teloclip.nf workflow

- `--bam`: Input BAM file (required)
- `--ref`: Reference FASTA file (required)
- `--outdir`: Output directory (default: "results")
- `--motifs`: Telomeric motif to search for (default: "TTAGGG")
- `--min_repeats`: Minimum consecutive motif repeats (default: 1)

## Output Structure

```text
results/
├── ref/
│   └── reference.fa.fai
├── filtered/
│   ├── overhangs.bam
│   └── overhangs.bam.bai
├── extracted/
│   └── overhangs/
│       ├── contig1_left.fa
│       ├── contig1_right.fa
│       └── ...
└── extended/
    ├── teloclip_extended.fasta
    └── teloclip_extension_report.txt
```

## Integration with nf-core

These modules follow nf-core guidelines and can be integrated into nf-core pipelines:

1. Place modules in your pipeline's `modules/local/` directory
2. Import and use as shown above
3. Configure via `nextflow.config`
4. Test with nf-core tools

## Requirements

- Nextflow >= 21.04.0
- Docker or Singularity
- Input BAM file with aligned long reads
- Reference FASTA file

## Troubleshooting

### Container not found

Ensure Docker is running and you have internet access to pull images.

### Permission errors

Check file permissions and ensure Docker has access to input/output directories.

### Out of memory

Increase resources in your nextflow.config:

```groovy
process {
    withName: TELOCLIP_FILTER {
        memory = '8 GB'
    }
}
```

## References

- [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html)
- [nf-core modules](https://nf-co.re/modules)
- [Teloclip documentation](https://github.com/adamtaranto/teloclip)
