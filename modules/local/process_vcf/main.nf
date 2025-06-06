process PROCESS_VCF {
    label 'process_single'
    tag "${meta.id}"

    // Using artic as it has both bcftools and pysam required here
    //  May look into a different container later
    // process_gvcf.py is from https://github.com/jts/ncov2019-artic-nf/blob/master/bin/process_gvcf.py with some small updates
    //  Notably adding in the genotype format as 1 for variants to work with newer bcftools, a slight qual filter, and tsv output to table later
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.6.2--pyhdfd78af_0' :
        'biocontainers/artic:1.6.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("${meta.id}.consensus.norm.vcf.gz"), path("${meta.id}.consensus.norm.vcf.gz.tbi"), emit: consensus_vcf
    tuple val(meta), path("${meta.id}.ambiguous.norm.vcf.gz"), path("${meta.id}.ambiguous.norm.vcf.gz.tbi"), emit: ambiguous_vcf
    tuple val(meta), path("${meta.id}.variants.vcf"), emit: variants_vcf
    tuple val(meta), path("${meta.id}.consensus.tsv"), emit: variants_tsv
    path "versions.yml", emit: versions

    script:
    def frameshiftArg = ''
    if ( params.no_frameshifts ) {
        frameshiftArg = '--no-frameshifts'
    }
    """
    # Process the VCF
    process_vcf.py \\
        -d ${params.min_depth} \\
        -l ${params.min_ambiguity_threshold} \\
        -u ${params.max_ambiguity_threshold} \\
        -q 20 \\
        -c ${meta.id}.processed.vcf \\
        -v ${meta.id}.variants.vcf \\
        -t ${meta.id}.consensus.tsv \\
        $frameshiftArg \\
        $vcf

    # Normalize variant records into canonical VCF representation
    bcftools norm \\
        -f $reference \\
        ${meta.id}.processed.vcf \\
        > ${meta.id}.processed.norm.vcf

    for vt in "ambiguous" "consensus"; do
        cat ${meta.id}.processed.norm.vcf \\
            | awk -v vartag=ConsensusTag=\$vt '\$0 ~ /^#/ || \$0 ~ vartag' > ${meta.id}.\$vt.norm.vcf

        bgzip -f ${meta.id}.\$vt.norm.vcf
        tabix -p vcf ${meta.id}.\$vt.norm.vcf.gz
    done

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        x:y
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.consensus.norm.vcf.gz
    touch ${meta.id}.consensus.norm.vcf.gz.tbi

    touch ${meta.id}.ambiguous.norm.vcf.gz
    touch ${meta.id}.ambiguous.norm.vcf.gz.tbi

    touch ${meta.id}.variants.vcf

    touch ${meta.id}.consensus.tsv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        x:y
    END_VERSIONS
    """
}
