// Maybe there is a better way to adjust the header than additional val inputs
//  But this works pretty well
process ADJUST_FASTA_HEADER {
    label 'process_single'
    tag "$meta.id"

    container "biocontainers/coreutils:8.31--h14c3975_0"

    input:
    tuple val(meta), path(consensus)
    tuple val(meta2), path(reference)
    val additional_extension_str
    val additional_header_str

    output:
    tuple val(meta), path("${meta.id}*.fasta"), emit: consensus
    path "versions.yml", emit: versions

    script:
    """
    sed "1 s/.*/>${meta.id}${additional_header_str} ${meta2.id}/" $consensus > ${meta.id}${additional_extension_str}.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version) | head -n 1 | cut -d' ' -f 4)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.consensus.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version) | head -n 1 | cut -d' ' -f 4)
    END_VERSIONS
    """
}
