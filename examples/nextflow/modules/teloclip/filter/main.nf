process TELOCLIP_FILTER {
    tag "$meta.id"
    label 'process_medium'

    container 'adamtaranto/teloclip:latest'

    input:
    tuple val(meta), path(sam)
    path fai

    output:
    tuple val(meta), path("*_overhangs.sam"), emit: sam
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    teloclip filter \\
        --ref-idx ${fai} \\
        ${args} \\
        ${sam} > ${prefix}_overhangs.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        teloclip: \$(teloclip --version | sed 's/teloclip, version //')
    END_VERSIONS
    """
}
