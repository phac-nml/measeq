// Two processes, big difference is the rg taken out from custom one
process ARTIC_MAKE_DEPTH_MASK{
    label 'process_small'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.10.3--pyhdfd78af_0' :
        'biocontainers/artic:1.10.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("${meta.id}.coverage_mask.txt"), emit: coverage_mask
    path "versions.yml", emit: versions

    script:
    """
    artic_make_depth_mask \\
        --depth ${params.min_depth} \\
        --store-rg-depths \\
        $reference \\
        $bam \\
        ${meta.id}.coverage_mask.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.coverage_mask.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}

// Slow but the bedtools adaptation I was working on I couldn't quite get to be genomic index
//  Will have to look at that more as it was a lot quicker
process CUSTOM_MAKE_DEPTH_MASK {
    label 'process_small'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.10.3--pyhdfd78af_0' :
        'biocontainers/artic:1.10.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("${meta.id}.coverage_mask.txt"), emit: coverage_mask
    path "versions.yml", emit: versions

    script:
    """
    cs_make_depth_mask.py \\
        --depth ${params.min_depth} \\
        $reference \\
        $bam \\
        ${meta.id}.coverage_mask.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.coverage_mask.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
