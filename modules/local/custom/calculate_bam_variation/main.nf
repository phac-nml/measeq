process CALCULATE_BAM_VARIATION {
    // High memory as not currently multi-threaded and some files are big
    label 'process_high_memory'
    tag "$meta.id"

    // I'm just going to use artic for now as it has updated dependencies and is used elsewhere
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(reference)

    output:
    tuple val(meta), path("${meta.id}_variation.csv"), optional: true, emit: csv
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    calc_bam_variation.py \\
        --bam $bam \\
        --reference $reference \\
        $args \\
        --sample ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        calc_bam_variation.py: 0.1.0
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_variation.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        calc_bam_variation.py: 0.1.0
    END_VERSIONS
    """
}
