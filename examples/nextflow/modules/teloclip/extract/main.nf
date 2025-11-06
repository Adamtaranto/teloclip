process TELOCLIP_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    container 'adamtaranto/teloclip:latest'

    input:
    tuple val(meta), path(sam)
    path fai

    output:
    tuple val(meta), path("${prefix}/*"), emit: fastas
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    teloclip extract \\
        --ref-idx ${fai} \\
        --extract-dir ${prefix} \\
        ${args} \\
        ${sam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        teloclip: \$(teloclip --version | sed 's/teloclip, version //')
    END_VERSIONS
    """
}
