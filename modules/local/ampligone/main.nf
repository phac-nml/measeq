process AMPLIGONE {
    label 'process_medium'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampligone:2.0.1--pyhdfd78af_0' :
        'biocontainers/ampligone:2.0.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastqs)
    tuple val(meta2), path(reference)
    path primers

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "detected_primers.bed", emit: primeres, optional: true
    path "versions.yml", emit: versions

    script:
    if (meta.single_end) {
        """
        ampligone \\
            --input $fastqs \\
            --output ${meta.id}.primertrimmed.fastq \\
            --reference $reference \\
            --primers $primers \\
            --amplicon-type end-to-end \\
            --export-primers detected_primers.bed \\
            -to \\
            --threads $task.cpus

        gzip ${meta.id}.primertrimmed.fastq

        # Versions #
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            AmpliGone: X
        END_VERSIONS
        """
    } else {
        // Fix later, we don't want to just use .fastp but for testing
        """
        for fastq in $fastqs; do
            if echo "\$fastq" | grep -qi '_1.fastp.fastq.gz'; then
                outext="_1.fastq"
            else
                outext="_2.fastq"
            fi

            ampligone \\
                --input \$fastq \\
                --output ${meta.id}_ampligone\$outext \\
                --reference $reference \\
                --primers $primers \\
                --amplicon-type end-to-end \\
                --export-primers primers.bed \\
                --threads $task.cpus

            gzip ${meta.id}_ampligone\$outext
        done

        # Versions #
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            AmpliGone: X
        END_VERSIONS
        """
    }

    stub:
    """
    touch ${meta.id}.coverage_mask.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AmpliGone: X
    END_VERSIONS
    """
}
