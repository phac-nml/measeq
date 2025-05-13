process MAKE_CUSTOM_REPORT {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker.io/darianhole/measeq-report:latest"

    input:
    path overall_qc_csv
    path depth_tsvs
    path n_depth_tsvs
    path baseq_tsvs
    val strain
    path report_template
    path subpages

    output:
    path "*.html", emit: html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Organize a bit
    mkdir -p positional_depth
    mv $depth_tsvs positional_depth/

    mkdir -p positional_n_depth
    mv $n_depth_tsvs positional_n_depth/

    mkdir -p positional_baseq
    mv $baseq_tsvs positional_baseq/

    # Create Report #
    Rscript -e "rmarkdown::render('$report_template', params = list(strain = '$strain', overall_qc = '$overall_qc_csv'))"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        x: 1
    END_VERSIONS
    """

    stub:
    """
    touch MeaSeq_Report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        x: 1
    END_VERSIONS
    """
}
