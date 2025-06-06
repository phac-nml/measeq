process MAKE_CUSTOM_REPORT {
    label 'process_medium'

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
    val genotype
    path report_template
    path subpages

    output:
    path "*.html", emit: html
    path "*_files", emit: data_dir // Using self-contained false 
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

    mkdir -p variation
    mv $variation_csvs variation/

    mkdir -p variant_tsv
    mv $variants_tsv variant_tsv/

    # Create Report #
    Rscript -e "rmarkdown::render('$report_template', params = list(genotype = '$genotype', overall_qc = '$overall_qc_csv'))"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        x:y
    END_VERSIONS
    """

    stub:
    """
    touch MeaSeq_Report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        x:y
    END_VERSIONS
    """
}
