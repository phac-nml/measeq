process VCF_TO_TSV {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.8.5--pyhdfd78af_0' :
        'biocontainers/artic:1.8.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.consensus.tsv"), optional: true, emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    vcf_to_tsv.py \\
        --sample "${meta.id}" \\
        --vcf $vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
        vcf_to_tsv: 0.1.0
    END_VERSIONS
    """

    stub:
    """
    ${meta.id}.consensus.tsv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
        vcf_to_tsv: 0.1.0
    END_VERSIONS
    """
}
