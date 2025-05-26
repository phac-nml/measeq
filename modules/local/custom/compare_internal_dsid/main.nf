/*
Exact matching of the N450 sequence to input DSID fasta file to determine if we already have the ID
    User has to provide the DSID fasta organized as:
    >DSID
    GTCAGTTCCA...

    Output is just a TSV that has:
    sample   matched_dsid
*/
process COMPARE_INTERNAL_DSID {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.10.0--h9ee0642_0' :
        'biocontainers/seqkit:2.10.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(n450)
    path id_fasta

    output:
    tuple val(meta), path("${meta.id}_id.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    result=""

    # If there is data in the file we want to match it, otherwise we just output the default None
    if [ "\$(wc -l < $n450)" -gt 1 ]; then
        # Create pattern from input fasta file
        pattern=\$(tail -n +2 $n450 | tr -d '\\n')

        # Search for a match
        result=\$(seqkit grep -sp \$pattern -m 0 $id_fasta | head -n 1 | sed -n 's/^>\\([0-9][0-9]*\\).*/\\1/p')
    fi

    # Output
    echo -e "sample\\tmatched_dsid" > ${meta.id}_id.tsv
    echo -e "${meta.id}\\t\${result:-None}" >> ${meta.id}_id.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(echo \$(seqkit version | sed 's/^seqkit v//'))
    END_VERSIONS
    """

    stub:
    """
    ${meta.id}_id.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(echo \$(seqkit version | sed 's/^seqkit v//'))
    END_VERSIONS
    """
}
