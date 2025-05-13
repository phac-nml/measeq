process BCFTOOLS_NORM {
    label 'process_small'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5acacb55c52bec97c61fd34ffa8721fce82ce823005793592e2a80bf71632cd0/data':
        'community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11' }"

    input:
    tuple val(meta), path(preconsensus), path(pass_vcf), path(pass_vcf_tbi)

    output:
    tuple val(meta), path("${meta.id}.pass.norm.vcf.gz"), path("${meta.id}.pass.norm.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    # Fixes variants that are in both the pass and fail vcf that were masked #
    bcftools norm \\
        --check-ref x \\
        -f $preconsensus \\
        $pass_vcf \\
        > ${meta.id}.pass.norm.vcf

    bgzip ${meta.id}.pass.norm.vcf
    tabix ${meta.id}.pass.norm.vcf.gz

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.pass.norm.vcf.gz
    touch ${meta.id}.pass.norm.vcf.gz.tbi

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
