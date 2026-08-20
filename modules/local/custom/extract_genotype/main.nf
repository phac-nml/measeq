process EXTRACT_GENOTYPE {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/coreutils:8.31--h14c3975_0'
        : 'biocontainers/coreutils:8.31--h14c3975_0' }"

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), env('GENOTYPE'), emit: genotype

    script:
    """
    col=\$(awk -F';' 'NR==1{for(i=1;i<=NF;i++){if(\$i=="clade"){print i; exit}}}' "$csv")
    GENOTYPE=\$(awk -F';' -v c="\$col" 'NR==2{print \$c; exit}' "$csv")
    """

    stub:
    """
    echo "default"
    """
}
