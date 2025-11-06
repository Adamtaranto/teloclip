process TELOCLIP_EXTEND {
    tag "$meta.id"
    label 'process_medium'

    container 'adamtaranto/teloclip:latest'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}_extended.fasta"), emit: fasta
    tuple val(meta), path("${prefix}_extension_report.txt"), emit: report
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    teloclip extend \\
        ${bam} \\
        ${fasta} \\
        --output-fasta ${prefix}_extended.fasta \\
        --stats-report ${prefix}_extension_report.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        teloclip: \$(teloclip --version | sed 's/teloclip, version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_extended.fasta
    touch ${prefix}_extension_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        teloclip: \$(teloclip --version | sed 's/teloclip, version //')
    END_VERSIONS
    """
}
