//
// Subworkflow for amplicon and non-amplicon consensus sequence generation for ONT data
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Initial Steps
include { ARTIC_GET_MODELS        } from '../../../modules/local/artic/get_model/main'
include { NANOQ                   } from '../../../modules/nf-core/nanoq/main'
include { MINIMAP2_ALIGN          } from '../../../modules/local/minimap2/main'
// Amplicon Variant Calling
include { ARTIC_ALIGN_TRIM        } from '../../../modules/local/artic/subcommands/main'
include { CLAIR3_POOL             } from '../../../modules/local/clair3/main'
include { ARTIC_VCF_MERGE         } from '../../../modules/local/artic/subcommands/main'
include { ARTIC_MAKE_DEPTH_MASK   } from '../../../modules/local/artic/subcommands/main'
// Non-Amplicon Variant Calling
include { CLAIR3_NO_POOL          } from '../../../modules/local/clair3/main'
include { CUSTOM_MAKE_DEPTH_MASK  } from '../../../modules/local/artic/subcommands/main'
// Both Post Variant Calling
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

workflow NANOPORE_CONSENSUS {

    take:
    ch_reference           // channel: [ [id], fasta ]
    ch_reference_fai       // channel: [ fai ]
    ch_input_fastqs        // channel: [ [id, single_end], fastqs ]
    ch_primer_bed          // channel: [ bed ]
    ch_split_amp_pools_bed // channel: [ pool, bed ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Model download if we are not using a local one
    //
    if ( params.local_model ) {
        ch_model = Channel.value(file(params.local_model, type: 'dir', checkIfExists: true))
    } else {
        ARTIC_GET_MODELS(
            params.model
        )
        ch_model = ARTIC_GET_MODELS.out.model
        ch_versions = ch_versions.mix(ARTIC_GET_MODELS.out.versions)
    }

    //
    // MODULE: Run Nanoq for input fastq quality filtering
    //
    NANOQ(
        ch_input_fastqs,
        'fastq'
    )
    ch_versions = ch_versions.mix(NANOQ.out.versions.first())

    //
    // MODULE: Run Minimap2
    //
    MINIMAP2_ALIGN(
        NANOQ.out.reads,
        ch_reference
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    //
    // PROCESS: Amplicon specific steps
    //
    if ( ch_primer_bed ) {
        //
        // MODULE: Align trim
        //  Trimming based on primers as there are a lot and don't want them
        //  affecting the variants
        ARTIC_ALIGN_TRIM (
            MINIMAP2_ALIGN.out.bam_bai,
            ch_primer_bed,
            'primers'
        )
        ch_versions = ch_versions.mix(ARTIC_ALIGN_TRIM.out.versions.first())
        ch_trimmed_bams_w_pool = ARTIC_ALIGN_TRIM.out.bam.combine(ch_split_amp_pools_bed)

        //
        // MODULE: Clair3 with pools
        //  Each pools run separately for variant calls 
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

        //
        // MODULE: Merge pooled VCFs
        //  To merge vcfs, have to utilize the transformVCFList function based on how artic handles input
        ARTIC_VCF_MERGE(
            ch_pooled_vcfs,
            ch_primer_bed
        )
        ch_versions = ch_versions.mix(ARTIC_VCF_MERGE.out.versions.first())
        ch_vcf = ARTIC_VCF_MERGE.out.vcf

        //
        // MODULE: Make depth mask based on minimum depth to call position
        //
        ARTIC_MAKE_DEPTH_MASK(
            ARTIC_ALIGN_TRIM.out.bam,
            ch_reference
        )
        ch_versions = ch_versions.mix(ARTIC_MAKE_DEPTH_MASK.out.versions.first())
        ch_depth_mask = ARTIC_MAKE_DEPTH_MASK.out.coverage_mask
    } else {
        //
        // PROCESS: Non-Amplicon based data
        //

        //
        // MODULE: Clair3 with no pool splitting
        //
        CLAIR3_NO_POOL(
            MINIMAP2_ALIGN.out.bam_bai,
            ch_reference,
            ch_reference_fai,
            ch_model
        )
        ch_versions = ch_versions.mix(CLAIR3_NO_POOL.out.versions.first())
        ch_vcf = CLAIR3_NO_POOL.out.vcf

        //
        // MODULE: Make depth mask based on minimum depth to call position
        //
        CUSTOM_MAKE_DEPTH_MASK(
            MINIMAP2_ALIGN.out.bam_bai,
            ch_reference
        )
        ch_versions = ch_versions.mix(CUSTOM_MAKE_DEPTH_MASK.out.versions.first())
        ch_depth_mask = CUSTOM_MAKE_DEPTH_MASK.out.coverage_mask

    }

    //
    // MODULE: Filter VCF based on input parameters
    //
    CUSTOM_VCF_FILTER(
        ch_vcf
    )
    ch_versions = ch_versions.mix(CUSTOM_VCF_FILTER.out.versions.first())

    //
    // MODULE: Masl failing regions
    //
    ARTIC_MASK(
        ch_depth_mask.join(CUSTOM_VCF_FILTER.out.fail_vcf, by: [0]),
        ch_reference
    )
    ch_versions = ch_versions.mix(ARTIC_MASK.out.versions.first())

    //
    // MODULE: Normalize
    //
    BCFTOOLS_NORM(
        ARTIC_MASK.out.preconsensus
            .join(CUSTOM_VCF_FILTER.out.pass_vcf, by: [0])
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

    //
    // MODULE: Make the final consensus
    //
    BCFTOOLS_CONSENSUS(
        BCFTOOLS_NORM.out.vcf
            .join(ARTIC_MASK.out.preconsensus, by: [0])
            .join(ch_depth_mask, by: [0])
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first())

    //
    // MODULE: Adjust the final fasta header for easier downstream analysis
    //
    ADJUST_FASTA_HEADER(
        BCFTOOLS_CONSENSUS.out.fasta,
        ch_reference,
        '.consensus',
        ''
    )
    ch_versions = ch_versions.mix(ADJUST_FASTA_HEADER.out.versions.first())

    emit:
    nanoq_stats = NANOQ.out.stats
    consensus   = ADJUST_FASTA_HEADER.out.consensus
    bam_bai     = ARTIC_ALIGN_TRIM.out.bam
    vcf         = BCFTOOLS_NORM.out.vcf
    versions    = ch_versions
}
