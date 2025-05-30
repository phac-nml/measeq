process SAMPLE_AMPLICON_DEPTH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.30.0--h9ee0642_0' :
        'biocontainers/csvtk:0.30.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*_amplicon_depth.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    csvtk cut \\
        -tT \\
        -f 4,7 \\
        $bed \\
    | csvtk replace \\
        -tT \\
        -f 2 \\
        -p "^0\$" \\
        -r 0.1 \\
        > ${meta.id}_amplicon_depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_amplicon_depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
