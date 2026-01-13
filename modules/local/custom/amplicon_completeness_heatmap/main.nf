process AMPLICON_COMPLETENESS_HEATMAP {
    label 'process_single'
    tag "$ref_id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.30.0--h9ee0642_0' :
        'biocontainers/csvtk:0.30.0--h9ee0642_0' }"

    input:
    tuple val(ref_id), path(tsvs)

    output:
    tuple val(ref_id), path("${ref_id}_amplicon_completeness_heatmap_mqc.tsv"), emit: heatmap_tsv
    path "versions.yml", emit: versions

    script:
    // Could likely be a pipeline operator but for now using csvtk
    """
    csvtk concat -tT $tsvs > concat.tsv
    csvtk sort -tTk 1:N concat.tsv > ${ref_id}_amplicon_completeness_heatmap_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """

    stub:
    """
    touch ${ref_id}_amplicon_completeness_heatmap_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
