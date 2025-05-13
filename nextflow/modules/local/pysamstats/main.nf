process PYSAMSTATS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysamstats:1.1.2--py39he47c912_12':
        'biocontainers/pysamstats:1.1.2--py39he47c912_12' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val type // Values: https://github.com/alimanfoo/pysamstats?tab=readme-ov-file#usage

    output:
    tuple val(meta), path("${meta.id}.${type}.stats.tsv"), emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pysamstats \\
        --type $type \\
        -d $bam \\
        > ${meta.id}.${type}.stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysamstats: \$(pysamstats -h | tail -n 2 | grep -Eo ": \\S+" | cut -d" " -f2)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.${type}.stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysamstats: \$(pysamstats -h | tail -n 2 | grep -Eo ": \\S+" | cut -d" " -f2)
    END_VERSIONS
    """
}
