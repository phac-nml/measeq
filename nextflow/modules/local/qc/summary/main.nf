// Summary QC
//  Using artic env for the moment as it has a bunch of tools and
//  is used in earlier steps
process MAKE_FINAL_QC_CSV {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.6.2--pyhdfd78af_0' :
        'biocontainers/artic:1.6.2--pyhdfd78af_0' }"

    input:
    path concat_csv
    path metadata
    val neg_control_pct_threshold
    val neg_ctrl_substrings

    output:
    path "overall.qc.csv", emit: csv
    path "versions.yml", emit: versions

    script:
    def version = workflow.manifest.version
    def metadataArg = metadata ? "--metadata $metadata" : ""
    """
    summary_qc.py \\
        --csv $concat_csv \\
        $metadataArg \\
        --threshold $neg_control_pct_threshold \\
        --neg_ctrl_substrings '$neg_ctrl_substrings' \\
        --version $version

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch overall.qc.csv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
