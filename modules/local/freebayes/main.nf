process FREEBAYES {
    label 'process_medium'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.9--hbefcdb2_1' :
        'biocontainers/freebayes:1.3.9--hbefcdb2_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("${meta.id}.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    freebayes \\
        -b $bam \\
        -f $reference \\
        -p 1 \\
        -C 1 \\
        -F 0.05 \\
        --pooled-continuous \\
        --min-coverage 5 \\
        --standard-filters \\
        | sed s/QR,Number=1,Type=Integer/QR,Number=1,Type=Float/ > ${meta.id}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        x:x
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        x:x
    END_VERSIONS
    """
}
