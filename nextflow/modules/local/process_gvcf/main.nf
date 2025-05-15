process PROCESS_GVCF {
    label 'process_single'
    tag "${meta.id}"

    // Using artic as it has both bcftools and pysam required here
    //  May look into a different container later
    //  And potentially update the process_gvcf script as for getting things working it was used as is
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.6.2--pyhdfd78af_0' :
        'biocontainers/artic:1.6.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(gvcf)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("${meta.id}.fixed.norm.vcf.gz"), path("${meta.id}.fixed.norm.vcf.gz.tbi"), emit: consensus_vcf
    tuple val(meta), path("${meta.id}.ambiguous.norm.vcf.gz"), path("${meta.id}.ambiguous.norm.vcf.gz.tbi"), emit: ambiguous_vcf
    tuple val(meta), path("${meta.id}.variants.vcf"), emit: variants_vcf
    tuple val(meta), path("${meta.id}.mask.txt"), emit: mask
    path "versions.yml", emit: versions

    script:
    """
    # Process the GVCF
    process_gvcf.py \\
        -d 10 \\
        -l 0.25 \\
        -u 0.75 \\
        -m ${meta.id}.mask.txt \\
        -c ${meta.id}.consensus.vcf \\
        -v ${meta.id}.variants.vcf \\
        $gvcf

    # Normalize variant records into canonical VCF representation
    bcftools norm \\
        -f $reference \\
        ${meta.id}.consensus.vcf \\
        > ${meta.id}.consensus.norm.vcf 

    for vt in "ambiguous" "fixed"; do
        cat ${meta.id}.consensus.norm.vcf \\
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
    touch ${meta.id}.fixed.norm.vcf.gz
    touch ${meta.id}.fixed.norm.vcf.gz.tbi

    touch ${meta.id}.ambiguous.norm.vcf.gz
    touch ${meta.id}.ambiguous.norm.vcf.gz.tbi

    touch ${meta.id}.variants.vcf
    touch ${meta.id}.mask.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        x:y
    END_VERSIONS
    """
}
