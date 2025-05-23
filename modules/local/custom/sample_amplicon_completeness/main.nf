process SAMPLE_AMPLICON_COMPLETENESS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(consensus)
    path amplicon_bed

    output:
    tuple val(meta), path("*_amplicon_completeness.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    calc_amp_completeness.py \\
        --sample ${meta.id} \\
        --consensus $consensus \\
        --amplicon_bed $amplicon_bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        x:x
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_amplicon_depth.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        x:x
    END_VERSIONS
    """
}
