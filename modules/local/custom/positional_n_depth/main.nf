process POSITIONAL_N_DEPTH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}.per_base_n.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    samtools mpileup \\
        -f $reference \\
        --count-orphans \\
        --no-BAQ \\
        --ignore-RG \\
        $bam \\
        > ${meta.id}_output.pileup

    echo -e "position\\tN_density" \
        > ${meta.id}.per_base_n.tsv
    awk '{count = gsub("N", "N"); print \$2 "\\t" count}' \\
        ${meta.id}_output.pileup \\
        >> ${meta.id}.per_base_n.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.per_base_n.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
