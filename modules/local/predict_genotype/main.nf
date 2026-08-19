process PREDICT_GENOTYPE {
    label 'process_low'
    tag "$meta.id"

    // I'm just going to use artic for now as it has updated dependencies and is used elsewhere
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.10.3--pyhdfd78af_0' :
        'biocontainers/artic:1.10.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastqs)
    path(mmi)

    output:
    tuple val(meta), path(fastqs), env('GENOTYPE'), emit: samplesheet
    path "*.csv", emit: csv
    path "versions.yml", emit: versions

    script:
    def read1 = fastqs[0]
    def read2 = meta.single_end ? '' : "-2 ${fastqs[1]}"
    """
    predict_genotype.py -1 ${read1} ${read2} -s ${meta.id} -r ${mmi}

    GENOTYPE=\$(awk -F',' 'NR==2 {print \$2}' predictions.csv)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        predict_genotype.py: 0.1.0
    END_VERSIONS
    """

    stub:
    """
    echo "default"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        predict_genotype.py: 0.1.0
    END_VERSIONS
    """
}
