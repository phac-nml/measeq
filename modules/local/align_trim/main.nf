process ALIGN_TRIM {
    label 'process_medium'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    // No singularity container
    container "${ task.ext.override_configured_container_registry != false
        ? 'docker.io/ahmedabdalla4/align_trim-samtools:latest'
        : 'ahmedabdalla4/align_trim-samtools:latest' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path primer_bed
    val mode // 'primers' or 'start'

    output:
    tuple val(meta), path("${meta.id}.*trimmed.rg.sorted.bam"), path("${meta.id}.*trimmed.rg.sorted.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.amplicon_depths.tsv"), emit: amp_report
    tuple val(meta), path("${meta.id}.alignreport-${mode}.csv"), emit: align_report
    path "versions.yml", emit: versions

    script:
    def argsList = []
    if ( params.normalise_ont ) {
        argsList.add("--normalise ${params.normalise_ont}")
    }
    outName = "${meta.id}.trimmed.rg.sorted.bam"
    // Start mode = Trim to start of primers instead of ends
    if ( mode == "primers" ) {
        outName = "${meta.id}.primertrimmed.rg.sorted.bam"
    } else {
        argsList.add("--no-trim-primers")
    }
    // Remove or keep incorrect primers
    if ( params.ont_keep_incorrect_primers ) {
        argsList.add("--allow-incorrect-pairs")
    }
    def argsConfig = argsList.join(" ")
    """
    awk -F'\t' 'BEGIN{OFS=FS} {if (NF==6) print \$0, "ACCT"; else print \$0}' "${primer_bed}" > tmp.bed && mv tmp.bed "${primer_bed}"

    align_trim \\
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
        align_trim: \$(echo \$(align_trim --version 2>&1) | sed 's/align_trim //')
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
        align_trim: \$(echo \$(align_trim --version 2>&1) | sed 's/align_trim //')
    END_VERSIONS
    """
}
