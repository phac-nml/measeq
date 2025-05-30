// Individual QC
//  Using artic env for the moment as it has a bunch of tools and
//  is used in earlier steps
process MAKE_SAMPLE_QC_CSV {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.6.2--pyhdfd78af_0' :
        'biocontainers/artic:1.6.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(consensus), path(depth_bed), path(nextclade_n450), path(nextclade_full),
            path(vcf), path(tbi), path(read_json), path(dsid)
    val genotype
    path primer_bed

    output:
    tuple val(meta), path ("${meta.id}.qc.csv"), emit: csv
    path "versions.yml", emit: versions

    script:
    // Add in which read filtering program was used based on the name
    def readJsonArg = ""
    if (read_json.name.contains(".fastp.")) {
        readJsonArg = "--fastp_json $read_json"
    } else {
        readJsonArg = "--nanoq_json $read_json"
    }
    // Add matched dsid if we have it
    def dsidArg = dsid ? "--matched_dsid $dsid" : ""
    // Add primer bed if we have it
    def seqPrimerArg = primer_bed ? "--seq_bed $primer_bed" : ""
    """
    sample_qc.py \\
        --bam $bam \\
        --consensus $consensus \\
        --depth $depth_bed \\
        --nextclade_n450 $nextclade_n450 \\
        --nextclade_custom $nextclade_full \\
        --genotype $genotype \\
        --vcf $vcf \\
        $readJsonArg \\
        $dsidArg \\
        $seqPrimerArg \\
        --sample $meta.id

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.qc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
