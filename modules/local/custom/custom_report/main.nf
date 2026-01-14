process MAKE_CUSTOM_REPORT {
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "docker.io/darianhole/measeq-report:latest"

    input:
    path overall_qc_csv
    path depth_tsvs
    path n_depth_tsvs
    path baseq_tsvs
    path variation_csvs
    path variants_tsv
    path normalized_depth_csv
    path report_template
    path subpages
    path version_yml
    val pipeline_version
    val revision
    val nf_version
    val name
    val email
    val phone
    val website

    output:
    path "*.html", emit: html
    path version_yml, includeInputs: true, emit: full_versions // So iridanext plugin can find it

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

    mkdir -p variation
    mv $variation_csvs variation/

    mkdir -p variant_tsv
    mv $variants_tsv variant_tsv/

    # Create Report #
    Rscript -e "rmarkdown::render(
        '$report_template',
        params = list(
            overall_qc = '$overall_qc_csv',
            version = '$pipeline_version',
            revision = '$revision',
            nf_version = '$nf_version',
            name = '$name',
            email = '$email',
            phone = '$phone',
            website = '$website'
        ))"
    """

    stub:
    """
    touch MeaSeq_Report.html
    """
}
