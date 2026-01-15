/*
    Artic subcommands and some custom adaptations of them to work with clair3 without
        explicitly having to use the minion.py workflow
        Needed for some of the user arguments and for non-amp based data
*/

// Helper function for combine VCFs in the format needed for artic merge
//  Files are always NAME.POOL.vcf based on previous process
//  So if the name has a . the size() - 2 hopefully gets the right number
def transformVCFList (inputList) {
    def transformedOutput = inputList.collect { vcf ->
        def name_split = vcf.name.split(/\./)
        def pool = name_split[name_split.size()-2]
        "${pool}:${vcf}"
    }.join(" ")
    return transformedOutput
}

// Subcommands start here
process ARTIC_ALIGN_TRIM {
    label 'process_medium'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path primer_bed
    val mode // 'primers' or 'start'

    output:
    tuple val(meta), path("${meta.id}.*trimmed.rg.sorted.bam"), path("${meta.id}.*trimmed.rg.sorted.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.amplicon_depths.tsv"), emit: amp_report
    path "versions.yml", emit: versions

    script:
    def argsList = []
    if ( params.normalise_ont ) {
        argsList.add("--normalise ${params.normalise_ont}")
    }
    outName = "${meta.id}.trimmed.rg.sorted.bam"
    // Start mode = Trim to start of primers instead of ends
    if ( mode == "primers" ) {
        outName = "${meta.id}.primertrimmed.rg.sorted.bam"
        argsList.add("--trim-primers")
    }
    def argsConfig = argsList.join(" ")
    """
    align_trim \\
        $argsConfig \\
        --remove-incorrect-pairs \\
        --report ${meta.id}.alignreport-${mode}.csv \\
        --amp-depth-report ${meta.id}.amplicon_depths.tsv \\
        --primer-match-threshold 15 \\
        $primer_bed \\
        < $bam \\
    | samtools sort -T ${meta.id} - -o $outName

    samtools index $outName

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.trimmed.rg.sorted.bam
    touch ${meta.id}.trimmed.rg.sorted.bam.bai
    touch ${meta.id}.amplicon_depths.tsv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}

process ARTIC_VCF_MERGE {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    // The vcf_tuples input is [[ path(vcf), val(pool) ], [...]]
    //   The path(vcf) is turned into a string of the full path using the val() input type
    //   The process still works, just is a bit iffy I'd say
    input:
    tuple val(meta), path(vcfs)
    path primer_bed

    output:
    tuple val(meta), path("${meta.id}.merged.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def vcfs_in_str = transformVCFList(vcfs)
    """
    artic_vcf_merge \\
        ${meta.id} \\
        $primer_bed \\
        2> ${meta.id}.primersitereport.txt \\
        $vcfs_in_str

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.merged.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}

// Literally only custom as the real one doesn't have a quality argument
process CUSTOM_VCF_FILTER {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.pass.vcf.gz"), path("${meta.id}.pass.vcf.gz.tbi"), emit: pass_vcf
    tuple val(meta), path("${meta.id}.fail.vcf"), emit: fail_vcf
    path "versions.yml", emit: versions

    script:
    def frameshiftArg = ''
    if ( params.no_frameshifts ) {
        frameshiftArg = '--no-frameshifts'
    }
    """
    cs_vcf_filter.py \\
        $frameshiftArg \\
        --min-depth ${params.min_depth} \\
        --min-qual ${params.min_variant_qual_c3} \\
        $vcf \\
        ${meta.id}.pass.vcf \\
        ${meta.id}.fail.vcf
    bgzip -f ${meta.id}.pass.vcf
    tabix -p vcf ${meta.id}.pass.vcf.gz

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.pass.vcf.gz
    touch ${meta.id}.pass.vcf.gz.tbi
    touch ${meta.id}.fail.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}

process ARTIC_MAKE_DEPTH_MASK{
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("${meta.id}.coverage_mask.txt"), emit: coverage_mask
    path "versions.yml", emit: versions

    script:
    """
    artic_make_depth_mask \\
        --depth ${params.min_depth} \\
        --store-rg-depths \\
        $reference \\
        $bam \\
        ${meta.id}.coverage_mask.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.coverage_mask.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
// Slow but the bedtools adaptation I was working on I couldn't quite get to be genomic index
//  Will have to look at that more as it was a lot quicker
process CUSTOM_MAKE_DEPTH_MASK {
    label 'process_single'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("${meta.id}.coverage_mask.txt"), emit: coverage_mask
    path "versions.yml", emit: versions

    script:
    """
    cs_make_depth_mask.py \\
        --depth ${params.min_depth} \\
        $reference \\
        $bam \\
        ${meta.id}.coverage_mask.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.coverage_mask.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}

process ARTIC_MASK {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(coverage_mask), path(fail_vcf)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("${meta.id}.preconsensus.fasta"), emit: preconsensus
    path "versions.yml", emit: versions

    script:
    """
    artic_mask \\
        $reference \\
        $coverage_mask \\
        $fail_vcf \\
        ${meta.id}.preconsensus.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.preconsensus.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}

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
