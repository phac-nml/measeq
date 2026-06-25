//
// Two separate Clair3 processes for now, probably will combine later into 1
//  Separated as the one takes in the Pool Bed file while the other doesn't
//  To combine later, just have to add a check for the input probably
//
process CLAIR3_POOL {
    label 'process_medium'
    label 'error_retry'
    tag "${meta.id}-Pool${pool}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:2.0.1--py311hbc58adc_0' :
        'biocontainers/clair3:2.0.1--py311hbc58adc_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(pool), path(pool_bed)
    tuple val(meta2), path(reference)
    path fai
    path model

    output:
    tuple val(meta), path("${meta.id}.${pool}.vcf"), val(pool), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    # Filter bam to be the pool only
    #  Otherwise the full alignment mode may see the differences in depth at the overlap
    #  region and not call specific SNPs
    samtools view -b -r $pool $bam -o ${pool}.sorted.bam
    samtools index ${pool}.sorted.bam

    # Run clair3 after
    run_clair3.sh \\
        --bam_fn="${pool}.sorted.bam" \\
        --bed_fn=$pool_bed \\
        --ref_fn="\$PWD/$reference" \\
        --threads=${task.cpus} \\
        --platform='ont' \\
        --model_path="$model" \\
        --output="${meta.id}-out" \\
        --haploid_sensitive \\
        --enable_long_indel \\
        --include_all_ctgs \\
        --chunk_size=5000 \\
        --ref_pct_full=1 \\
        --var_pct_full=1 \\
        --no_phasing_for_fa \\
        --enable_variant_calling_at_sequence_head_and_tail

    gunzip ${meta.id}-out/merge_output.vcf.gz
    ln -s ${meta.id}-out/merge_output.vcf ${meta.id}.${pool}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(run_clair3.sh -v) | sed 's/Clair3 //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.${pool}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(run_clair3.sh -v) | sed 's/Clair3 //')
    END_VERSIONS
    """
}

process CLAIR3_NO_POOL {
    label 'process_medium'
    label 'error_retry'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:2.0.1--py311hbc58adc_0' :
        'biocontainers/clair3:2.0.1--py311hbc58adc_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(reference)
    path fai
    path model

    output:
    tuple val(meta), path("${meta.id}.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    run_clair3.sh \\
        --bam_fn=$bam \\
        --ref_fn="\$PWD/$reference" \\
        --threads=${task.cpus} \\
        --platform='ont' \\
        --model_path="$model" \\
        --output="${meta.id}-out" \\
        --haploid_sensitive \\
        --enable_long_indel \\
        --include_all_ctgs \\
        --chunk_size=5000 \\
        --no_phasing_for_fa \\
        --ref_pct_full=1 \\
        --var_pct_full=1 \\
        --enable_variant_calling_at_sequence_head_and_tail

    gunzip ${meta.id}-out/merge_output.vcf.gz
    ln -s ${meta.id}-out/merge_output.vcf ${meta.id}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(run_clair3.sh -v) | sed 's/Clair3 //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(run_clair3.sh -v) | sed 's/Clair3 //')
    END_VERSIONS
    """
}
