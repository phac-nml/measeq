process ARTIC_MASK {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.10.3--pyhdfd78af_0' :
        'biocontainers/artic:1.10.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(coverage_mask), path(fail_vcf)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("${meta.id}.preconsensus.fasta"), emit: preconsensus
    path "versions.yml", emit: versions

    script:
    """
    artic_mask \\
        $reference \\
        $coverage_mask \\
        $fail_vcf \\
        ${meta.id}.preconsensus.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.preconsensus.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
