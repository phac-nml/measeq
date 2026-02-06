process GENERATE_REF_INTERMEDIATES {
    label 'process_single'
    tag "${meta.id}"

    conda "bioconda::samtools=1.19.2 bioconda::htslib=1.19.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0':
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(reference)

    output:
    tuple val(meta), path("${reference}.fai"), emit: fai
    tuple val(meta), path("*refstats.txt"), emit: refstats
    tuple val(meta), path("*genome.bed"), emit: genome_bed
    path "versions.yml", emit: versions

    script:
    """
    samtools faidx $reference
    cat ${reference}.fai | awk '{print \$1 ":1-" \$2+1}' > ${meta.id}.refstats.txt
    cat ${reference}.fai | awk '{ print \$1 "	0	" \$2 }' > ${meta.id}.genome.bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${reference}.fai
    touch refstats.txt
    touch genome.bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process GENERATE_AMPLICON_BED {
    label 'process_single'
    tag "${meta.id}"

    conda "conda-forge::python=3.10.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python%:3.10.2'
        : 'biocontainers/python:3.10.2' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*amplicon.bed"), emit: bed
    path "versions.yml", emit: versions

    script:
    """
    primers_to_amplicons.py \\
        --bed $bed

    mv amplicon.bed ${meta.id}.amplicon.bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        primers_to_amplicons: 0.1.0
    END_VERSIONS
    """

    stub:
    """
    touch amplicon.bed
    touch tiling_region.bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        primers_to_amplicons: 0.1.0
    END_VERSIONS
    """
}

process SPLIT_AMPLICON_REGION {
    label 'process_single'
    tag "${meta.id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/coreutils:8.31--h14c3975_0'
        : 'biocontainers/coreutils:8.31--h14c3975_0' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed"), emit: bed

    script:
    """
    awk -F'\t' -v OFS='\t' 'NR>0{print \$1, \$2, \$3, \$4, \$5, \$6 > \$5".${meta.id}.bed"}' $bed
    """

    stub:
    """
    touch 1.bed
    touch 2.bed
    """
}
