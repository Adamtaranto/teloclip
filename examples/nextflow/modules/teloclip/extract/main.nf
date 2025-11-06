process TELOCLIP_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    container 'adamtaranto/teloclip:latest'

    input:
    tuple val(meta), path(bam)
    path fai

    output:
    tuple val(meta), path("overhangs/*"), emit: fastas
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p overhangs

    teloclip extract \\
        --ref-idx ${fai} \\
        --outdir overhangs \\
        --prefix ${prefix} \\
        ${args} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        teloclip: \$(teloclip --version | sed 's/teloclip //')
    END_VERSIONS
    """
}
