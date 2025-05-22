process MINIMAP2_ALIGN {
    label 'process_medium'
    tag "$meta.id"

    // Using the seqera one as it has more updated tools that the mulled biocontainer one from what I can see
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/66/66dc96eff11ab80dfd5c044e9b3425f52d818847b9c074794cf0c02bfa781661/data' :
        'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c' }"

    input:
    tuple val(meta), path(fastq)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam_bai
    path "versions.yml", emit: versions

    script:
    """
    minimap2 \\
        -a \\
        -x map-ont \\
        -t task.cpus \\
        $reference \\
        $fastq \\
    | samtools view -b -F 4 - \\
    | samtools sort -o ${meta.id}.sorted.bam

    samtools index ${meta.id}.sorted.bam

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
