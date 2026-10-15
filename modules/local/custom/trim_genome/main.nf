process TRIM_GENOME_UTRS {
    label 'process_single'
    tag "$meta.id"

    // I'm just going to use artic for now as it has updated dependencies and is used elsewhere
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.10.3--pyhdfd78af_0' :
        'biocontainers/artic:1.10.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(trim_guide_fasta)

    output:
    tuple val(meta), path("*.trimmed.consensus.fasta"), optional: true, emit: trimmed
    path "versions.yml", emit: versions

    script:
    """
    trim_genome.py \
        --trim_guide $trim_guide_fasta \
        --input $fasta \
        --outfile ${meta.id}.trimmed.consensus.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        python: \$(python --version | sed 's/Python //g')
        trim_genome.py: 0.1.0
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.trimmed.consensus.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        python: \$(python --version | sed 's/Python //g')
        trim_genome.py: 0.1.0
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
