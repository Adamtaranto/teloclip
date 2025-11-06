#!/usr/bin/env nextflow

/*
 * Teloclip Example Workflow
 *
 * This workflow demonstrates how to use teloclip in a Nextflow pipeline
 * with proper containerization (each tool in its own container).
 *
 * Example usage:
 *   nextflow run teloclip.nf --bam input.bam --ref ref.fa
 */

nextflow.enable.dsl=2

// Parameters
params.bam = null
params.ref = null
params.outdir = "results"
params.motifs = "TTAGGG"
params.min_repeats = 1

// Check required parameters
if (!params.bam) error "Please provide --bam parameter"
if (!params.ref) error "Please provide --ref parameter"

/*
 * Process: Index reference FASTA
 */
process INDEX_FASTA {
    container 'quay.io/biocontainers/samtools:latest'
    publishDir "${params.outdir}/ref", mode: 'copy'

    input:
    path ref

    output:
    path "${ref}.fai"

    script:
    """
    samtools faidx ${ref}
    """
}

/*
 * Process: Convert BAM to SAM
 */
process BAM_TO_SAM {
    container 'quay.io/biocontainers/samtools:latest'

    input:
    path bam

    output:
    path "input.sam"

    script:
    """
    samtools view -h ${bam} > input.sam
    """
}

/*
 * Process: Filter alignments with teloclip
 */
process TELOCLIP_FILTER {
    container 'adamtaranto/teloclip:latest'
    publishDir "${params.outdir}/filtered", mode: 'copy'

    input:
    path sam
    path fai

    output:
    path "overhangs.sam"

    script:
    def motif_arg = params.motifs ? "--motifs ${params.motifs}" : ""
    def repeats_arg = params.min_repeats > 1 ? "--min-repeats ${params.min_repeats}" : ""
    """
    teloclip filter \
        --ref-idx ${fai} \
        ${motif_arg} \
        ${repeats_arg} \
        ${sam} > overhangs.sam
    """
}

/*
 * Process: Convert SAM to sorted BAM
 */
process SAM_TO_BAM {
    container 'quay.io/biocontainers/samtools:latest'
    publishDir "${params.outdir}/filtered", mode: 'copy'

    input:
    path sam

    output:
    tuple path("overhangs.bam"), path("overhangs.bam.bai")
    path sam

    script:
    """
    samtools sort ${sam} > overhangs.bam
    samtools index overhangs.bam
    """
}

/*
 * Process: Extract overhang sequences (optional)
 */
process TELOCLIP_EXTRACT {
    container 'adamtaranto/teloclip:latest'
    publishDir "${params.outdir}/extracted", mode: 'copy'

    input:
    path sam
    path fai

    output:
    path "overhangs/*.fasta"

    script:
    """
    mkdir -p overhangs
    teloclip extract \
        --ref-idx ${fai} \
        --extract-dir overhangs \
        ${sam}
    """
}

/*
 * Process: Extend telomeric sequences
 */
process TELOCLIP_EXTEND {
    container 'adamtaranto/teloclip:latest'
    publishDir "${params.outdir}/extended", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    path ref
    path fai

    output:
    path "*_extended.fasta"
    path "*_extension_report.txt"

    script:
    """
    teloclip extend \\
        ${bam} \\
        ${ref} \\
        --output-fasta teloclip_extended.fasta \\
        --stats-report teloclip_extension_report.txt \\
        --count-motifs ${params.motifs} \\
        --screen-terminal-bases 1000
    """
}

/*
 * Main workflow
 */
workflow {
    // Create channels
    bam_ch = Channel.fromPath(params.bam)
    ref_ch = Channel.fromPath(params.ref)

    // Index reference
    fai_ch = INDEX_FASTA(ref_ch)

    // Convert BAM to SAM
    sam_ch = BAM_TO_SAM(bam_ch)

    // Filter with teloclip
    filtered_sam = TELOCLIP_FILTER(sam_ch, fai_ch)

    // Convert to sorted and indexed BAM
    sam_bam_outputs = SAM_TO_BAM(filtered_sam)
    bam_indexed = sam_bam_outputs[0]
    filtered_sam_copy = sam_bam_outputs[1]

    // Extract sequences (optional - runs in parallel with extend)
    TELOCLIP_EXTRACT(filtered_sam_copy, fai_ch)

    // Extend reference genome with telomeric sequences
    TELOCLIP_EXTEND(bam_indexed, ref_ch, fai_ch)
}

workflow.onComplete {
    println "Pipeline completed!"
    println "Results saved to: ${params.outdir}"
    println "  - Filtered BAM: ${params.outdir}/filtered/"
    println "  - Extracted sequences: ${params.outdir}/extracted/"
    println "  - Extended genome: ${params.outdir}/extended/"
}
