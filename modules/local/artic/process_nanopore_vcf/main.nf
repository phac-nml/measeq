/*
    Custom previously to add in command line args for users to adjust in the pipeline
        Slight changes from artic where we've added:
            - A total depth check for each variant
            - Minimum quality skip to avoid unneeded masking
            - Minimum allele frequency to avoid unneeded masking
*/
process PROCESS_NANOPORE_VCF {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.10.3--pyhdfd78af_0' :
        'biocontainers/artic:1.10.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.pass.vcf.gz"), path("${meta.id}.pass.vcf.gz.tbi"), emit: pass_vcf
    tuple val(meta), path("${meta.id}.fail.vcf"), emit: fail_vcf
    tuple val(meta), path("${meta.id}.suppress.txt"), emit: suppress_txt
    path "versions.yml", emit: versions

    script:
    def frameshiftArg = ''
    if ( params.no_frameshifts ) {
        frameshiftArg = '--no-frameshifts'
    }
    """
    process_nanopore_vcf.py \\
        $frameshiftArg \\
        --min-depth ${params.min_depth} \\
        --min-qual ${params.min_variant_qual_c3} \\
        --min-frameshift-qual ${params.min_frameshift_qual_c3} \\
        --min-allele-freq ${params.min_allele_freq_c3} \\
        --min-mask-freq ${params.min_mask_freq_c3} \\
        --min-threshold-depth ${params.min_site_threshold_c3} \\
        $vcf \\
        $bam \\
        ${meta.id}.pass.vcf \\
        ${meta.id}.fail.vcf \\
        > ${meta.id}.suppress.txt
    bgzip -f ${meta.id}.pass.vcf
    tabix -p vcf ${meta.id}.pass.vcf.gz

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
        process_nanopore_vcf: 0.1.0
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
