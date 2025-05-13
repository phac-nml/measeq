process GENERATE_REF_INTERMEDIATES {
    label 'process_single'

    conda "bioconda::samtools=1.19.2 bioconda::htslib=1.19.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0':
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(reference)

    output:
    path "${reference}.fai", emit: fai
    path "refstats.txt", emit: refstats
    path "genome.bed", emit: genome_bed
    path "versions.yml", emit: versions

    script:
    """
    samtools faidx $reference
    cat ${reference}.fai | awk '{print \$1 ":1-" \$2+1}' > refstats.txt
    cat ${reference}.fai | awk '{ print \$1 "	0	" \$2 }' > genome.bed

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

process CREATE_AMPLICON_BED {
    label 'process_single'

    conda "conda-forge::python=3.10.2"
    container "quay.io/biocontainers/python:3.10.2"

    input:
    path bed

    output:
    path "amplicon.bed", emit: bed
    path "versions.yml", emit: versions

    script:
    """
    primers_to_amplicons.py \\
        --bed $bed

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
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
    END_VERSIONS
    """
}

process SPLIT_AMPLICON_REGION {
    label 'process_single'

    input:
    path bed

    output:
    path "*.bed", emit: bed

    script:
    """
    awk -F'\t' -v OFS='\t' 'NR>0{print \$1, \$2, \$3, \$4, \$5, \$6 > \$5".bed"}' $bed
    """

    stub:
    """
    touch 1.bed
    touch 2.bed
    """
}
