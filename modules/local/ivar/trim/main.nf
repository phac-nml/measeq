// This is basically the IVAR module from NF-Core with added samtools steps
//  This is to slightly simplify this part of the workflow
//  And if we want to add in the primer_pairs args
process IVAR_TRIM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ivar:1.4.4--h077b44d_0' :
        'biocontainers/ivar:1.4.4--h077b44d_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("${meta.id}.primertrimmed.sorted.bam"), path("${meta.id}.primertrimmed.sorted.bam.bai"), emit: bam
    tuple val(meta), path('*.log'), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ivar trim \\
        -m 30 -q 20 -e \\
        -i $bam \\
        -b $bed \\
        -p ${meta.id}.primertrimmed \\
        > ${meta.id}.ivar.log

    samtools \\
        sort \\
        -o ${meta.id}.primertrimmed.sorted.bam \\
        ${meta.id}.primertrimmed.bam
    samtools index ${meta.id}.primertrimmed.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version | sed -n 's|iVar version \\(.*\\)|\\1|p')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.ivar.log
    touch ${meta.id}.primertrimmed.sorted.bam
    touch ${meta.id}.primertrimmed.sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version | sed -n 's|iVar version \\(.*\\)|\\1|p')
    END_VERSIONS
    """
}
