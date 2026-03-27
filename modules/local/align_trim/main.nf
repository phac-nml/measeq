process ALIGN_TRIM {
    label 'process_medium'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path primer_bed
    val mode // 'primers' or 'start'

    output:
    tuple val(meta), path("${meta.id}.*trimmed.rg.sorted.bam"), path("${meta.id}.*trimmed.rg.sorted.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.amplicon_depths.tsv"), emit: amp_report, optional: true // only made if you normalize on 1.7.4
    tuple val(meta), path("${meta.id}.alignreport-${mode}.csv"), emit: align_report
    path "versions.yml", emit: versions

    script:
    def argsList = []
    argsList.add("--normalise ${params.normalise_ont}")
    outName = "${meta.id}.trimmed.rg.sorted.bam"
    // Start mode = Trim to start of primers instead of ends
    if ( mode == "primers" ) {
        outName = "${meta.id}.primertrimmed.rg.sorted.bam"
        argsList.add("--trim-primers")
    }
    // Remove or keep incorrect primers
    if ( ! params.ont_keep_incorrect_primers ) {
        argsList.add("--remove-incorrect-pairs")
    }
    def argsConfig = argsList.join(" ")
    """
    awk -F'\t' 'BEGIN{OFS=FS} {if (NF==6) print \$0, "NA"; else print \$0}' "${primer_bed}" > tmp.bed && mv tmp.bed "${primer_bed}"

    align_trim_seeded.py \\
        $argsConfig \\
        --report ${meta.id}.alignreport-${mode}.csv \\
        --amp-depth-report ${meta.id}.amplicon_depths.tsv \\
        --primer-match-threshold 15 \\
        $primer_bed \\
        < $bam \\
    | samtools sort -T ${meta.id} - -o $outName

    samtools index $outName

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.trimmed.rg.sorted.bam
    touch ${meta.id}.trimmed.rg.sorted.bam.bai
    touch ${meta.id}.amplicon_depths.tsv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
