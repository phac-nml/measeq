//
// Subworkflow for amplicon and non-amplicon consensus sequence generation
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Amplicon
include { ARTIC_ALIGN_TRIM        } from '../../../modules/local/artic/subcommands/main'
include { CLAIR3_POOL             } from '../../../modules/local/clair3/main'
include { ARTIC_VCF_MERGE         } from '../../../modules/local/artic/subcommands/main'
include { ARTIC_MAKE_DEPTH_MASK   } from '../../../modules/local/artic/subcommands/main'
// Non-Amplicon
include { CLAIR3_NO_POOL          } from '../../../modules/local/clair3/main'
include { CUSTOM_MAKE_DEPTH_MASK  } from '../../../modules/local/artic/subcommands/main'
// Both
include { CUSTOM_VCF_FILTER       } from '../../../modules/local/artic/subcommands/main'
include { ARTIC_MASK              } from '../../../modules/local/artic/subcommands/main'
include { BCFTOOLS_NORM           } from '../../../modules/local/bcftools/norm/main'
include { BCFTOOLS_CONSENSUS      } from '../../../modules/nf-core/bcftools/consensus/main'
include { ADJUST_FASTA_HEADER     } from '../../../modules/local/artic/subcommands/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO SETUP REFERENCE DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENERATE_CONSENSUS {

    take:
    ch_reference           //
    ch_reference_fai       //
    ch_bam_bai             //
    ch_model               //
    ch_primer_bed          // value: path to reference file or null
    ch_split_amp_pools_bed //

    main:
    ch_versions = Channel.empty()

    //
    // Amplicon
    //
    if ( ch_primer_bed ) {
        // Trimming based on primers as there are a lot and don't want them
        //  affecting the variants
        ARTIC_ALIGN_TRIM (
            ch_bam_bai,
            ch_primer_bed,
            'primers'
        )
        ch_versions = ch_versions.mix(ARTIC_ALIGN_TRIM.out.versions.first())
        ch_trimmed_bams_w_pool = ARTIC_ALIGN_TRIM.out.bam.combine(ch_split_amp_pools_bed)

        // Clair3 with pools
        CLAIR3_POOL(
            ch_trimmed_bams_w_pool,
            ch_reference,
            ch_reference_fai,
            ch_model
        )
        ch_versions = ch_versions.mix(CLAIR3_POOL.out.versions.first())

        // Merge pools by merging the vcf files for each pool together
        CLAIR3_POOL.out.vcf
            .map { it -> tuple(it[0], tuple(it[1], it[2])) }
            .groupTuple()
            .set { ch_pooled_vcfs } // Channel: [ val(meta), [[path(vcf), val(pool)], [...]] ]

        // To merge vcfs, have to utilize the transformVCFList function based on how artic handles input
        //  For now anyway, until a better solution comes up
        ARTIC_VCF_MERGE(
            ch_pooled_vcfs,
            ch_primer_bed
        )
        ch_versions = ch_versions.mix(ARTIC_VCF_MERGE.out.versions.first())
        ch_vcf = ARTIC_VCF_MERGE.out.vcf

        ARTIC_MAKE_DEPTH_MASK(
            ARTIC_ALIGN_TRIM.out.bam,
            ch_reference
        )
        ch_versions = ch_versions.mix(ARTIC_MAKE_DEPTH_MASK.out.versions.first())
        ch_depth_mask = ARTIC_MAKE_DEPTH_MASK.out.coverage_mask
    } else {
        //
        // Non-Amplicon
        //
        CLAIR3_NO_POOL(
            ch_bam_bai,
            ch_reference,
            ch_reference_fai,
            ch_model
        )
        ch_versions = ch_versions.mix(CLAIR3_NO_POOL.out.versions.first())
        ch_vcf = CLAIR3_NO_POOL.out.vcf

        CUSTOM_MAKE_DEPTH_MASK(
            ch_bam_bai,
            ch_reference
        )
        ch_versions = ch_versions.mix(CUSTOM_MAKE_DEPTH_MASK.out.versions.first())
        ch_depth_mask = CUSTOM_MAKE_DEPTH_MASK.out.coverage_mask

    }

    CUSTOM_VCF_FILTER(
        ch_vcf
    )
    ch_versions = ch_versions.mix(CUSTOM_VCF_FILTER.out.versions.first())

    ARTIC_MASK(
        ch_depth_mask.join(CUSTOM_VCF_FILTER.out.fail_vcf, by: [0]),
        ch_reference
    )
    ch_versions = ch_versions.mix(ARTIC_MASK.out.versions.first())

    BCFTOOLS_NORM(
        ARTIC_MASK.out.preconsensus
            .join(CUSTOM_VCF_FILTER.out.pass_vcf, by: [0])
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

    BCFTOOLS_CONSENSUS(
        BCFTOOLS_NORM.out.vcf
            .join(ARTIC_MASK.out.preconsensus, by: [0])
            .join(ch_depth_mask, by: [0])
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first())

    ADJUST_FASTA_HEADER(
        BCFTOOLS_CONSENSUS.out.fasta,
        ch_reference,
        '.consensus',
        ''
    )
    ch_versions = ch_versions.mix(ADJUST_FASTA_HEADER.out.versions.first())

    emit:
    consensus = ADJUST_FASTA_HEADER.out.consensus
    bam_bai   = ARTIC_ALIGN_TRIM.out.bam
    vcf       = BCFTOOLS_NORM.out.vcf
    versions  = ch_versions
}
