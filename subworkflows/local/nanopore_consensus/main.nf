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
include { ALIGN_TRIM              } from '../../../modules/local/artic/align_trim/main'
include { CLAIR3_POOL             } from '../../../modules/local/clair3/main'
include { ARTIC_VCF_MERGE         } from '../../../modules/local/artic/vcf_merge/main'
include { ARTIC_MAKE_DEPTH_MASK   } from '../../../modules/local/artic/make_depth_mask/main'
// Non-Amplicon Variant Calling
include { CLAIR3_NO_POOL          } from '../../../modules/local/clair3/main'
include { CUSTOM_MAKE_DEPTH_MASK  } from '../../../modules/local/artic/make_depth_mask/main'
// Both Post Variant Calling
include { PROCESS_NANOPORE_VCF    } from '../../../modules/local/artic/process_nanopore_vcf/main'
include { ARTIC_MASK              } from '../../../modules/local/artic/mask/main'
include { BCFTOOLS_NORM           } from '../../../modules/local/bcftools/norm/main'
include { BCFTOOLS_CONSENSUS      } from '../../../modules/nf-core/bcftools/consensus/main'
include { ADJUST_FASTA_HEADER     } from '../../../modules/local/adjust_fasta_header/main'
// For Final Reporting Only
include { VCF_TO_TSV              } from '../../../modules/local/custom/vcf_to_tsv/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO SETUP REFERENCE DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOPORE_CONSENSUS {

    take:
    ch_reference            // channel: [ meta_ref, fasta ]
    ch_fai                  // channel: [ meta_ref, fai ]
    ch_samples              // channel: [ meta[id, single_end, ref_id], fastqs ]
    ch_primer_bed           // channel: [ meta_ref, bed ]
    ch_split_amp_pools_bed  // channel: [ meta, bed ]

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

    // Check that we have the clair3 .pt file required and aren't giving an older model
    ch_model.map { model ->
        def pt_file = file("${model}/*.pt")
        if (!pt_file.exists()) {
            error("Missing .pt file in input model folder. Check that you are giving the path to a Clair3 v2.0.0+ pytorch model and not an older model")
        }
    }

    //
    // MODULE: Run Nanoq for input fastq quality filtering
    //
    NANOQ(
        ch_samples,
        'fastq'
    )
    ch_versions = ch_versions.mix(NANOQ.out.versions.first())

    //
    // MODULE: Run Minimap2
    //
    // Prepare Input
    ch_minimap_input = NANOQ.out.reads
        .map { meta, reads -> tuple(meta.ref_id, meta, reads) }
        .combine(ch_reference.map { meta_ref, fasta -> tuple(meta_ref.id, meta_ref, fasta) }, by: 0)
        .multiMap { _ref_id, meta, reads, meta_ref, fasta ->
            reads: tuple(meta, reads)
            reference: tuple(meta_ref, fasta)
        }

    // Run Module
    MINIMAP2_ALIGN(
        ch_minimap_input.reads,
        ch_minimap_input.reference
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())
    ch_bam_bai = MINIMAP2_ALIGN.out.bam_bai

    //
    // Prepare CLAIR3 reference file before if statement
    //
    ch_ref_with_fai = ch_reference
        .map { meta_ref, fasta -> tuple(meta_ref.id, meta_ref, fasta) }
        .join(ch_fai.map { meta_ref, fai -> tuple(meta_ref.id, fai) }, by: [0])
        .map { ref_id, meta_ref, fasta, fai ->
                tuple(ref_id, meta_ref, fasta, fai)
        }

    //
    // PROCESS: Amplicon specific steps
    //
    if ( params.amplicon || params.primer_bed ) {
        //
        // MODULE: Align trim
        //  Trimming based on primers as there are a lot and we don't want them
        //  affecting the variants
        //
        // Prepare Input
        ch_artic_align_input = ch_bam_bai
            .map { meta, bam, bai -> tuple(meta.ref_id, meta, bam, bai) }
            .combine(ch_primer_bed.map { meta_ref, bed -> tuple(meta_ref.id, bed) }, by: 0)
            .multiMap { _ref_id, meta, bam, bai, bed ->
                bam_bai: tuple(meta, bam, bai)
                bed: bed
            }

        // Run Module
        ALIGN_TRIM (
            ch_artic_align_input.bam_bai,
            ch_artic_align_input.bed,
            'primers'
        )
        ch_versions = ch_versions.mix(ALIGN_TRIM.out.versions.first())
        ch_bam = ALIGN_TRIM.out.bam

        //
        // MODULE: Clair3 with pools
        //  Each pools run separately for variant calls
        //
        // Prepare Inputs
        ch_trimmed_bams_w_pool = ch_bam
            .map { meta, bam, bai -> tuple(meta.ref_id, meta, bam, bai) }
            .combine(ch_split_amp_pools_bed.map { meta_ref_pool, pool_bed -> tuple(meta_ref_pool.id, meta_ref_pool.pool, pool_bed) }, by: 0)
            .map { _ref_id, meta, bam, bai, pool, pool_bed ->
                tuple(meta, bam, bai, pool, pool_bed)
            }

        // Combine input with reference data
        ch_clair3_inputs = ch_trimmed_bams_w_pool
            .map { meta, bam, bai, pool, pool_bed -> tuple(meta.ref_id, meta, bam, bai, pool, pool_bed) }
            .combine(ch_ref_with_fai, by: 0)
            .multiMap { _ref_id, meta, bam, bai, pool, pool_bed, meta_ref, fasta, fai ->
                bam_pool: tuple(meta, bam, bai, pool, pool_bed)
                reference: tuple(meta_ref, fasta)
                fai: fai
            }

        // Run Module
        CLAIR3_POOL(
            ch_clair3_inputs.bam_pool,
            ch_clair3_inputs.reference,
            ch_clair3_inputs.fai,
            ch_model
        )
        ch_versions = ch_versions.mix(CLAIR3_POOL.out.versions.first())

        // Merge pools by merging the vcf files for each pool together
        CLAIR3_POOL.out.vcf
            .map { meta, vcf, _pool -> tuple(meta, vcf) }
            .groupTuple()
            .set { ch_pooled_vcfs }

        //
        // MODULE: Merge pooled VCFs
        //  To merge vcfs, have to utilize the transformVCFList function based on how artic handles input
        //
        // Prepare Input
        ch_artic_vcf_input = ch_pooled_vcfs
            .map { meta, vcf -> tuple(meta.ref_id, meta, vcf) }
            .combine(ch_primer_bed.map { meta_ref, bed -> tuple(meta_ref.id, bed) }, by: 0)
            .multiMap { _ref_id, meta, vcf, bed ->
                vcf: tuple(meta, vcf)
                bed: bed
            }

        // Run Module
        ARTIC_VCF_MERGE(
            ch_artic_vcf_input.vcf,
            ch_artic_vcf_input.bed
        )
        ch_versions = ch_versions.mix(ARTIC_VCF_MERGE.out.versions.first())
        ch_vcf = ARTIC_VCF_MERGE.out.vcf

        //
        // MODULE: Make depth mask based on minimum depth to call position
        //
        // Prepare Inputs
        ch_artic_depth_input = ch_bam
            .map { meta, bam, bai -> tuple(meta.ref_id, meta, bam, bai) }
            .combine(ch_reference.map { meta_ref, fasta -> tuple(meta_ref.id, meta_ref, fasta) }, by: 0)
            .multiMap { _ref_id, meta, bam, bai, meta_ref, fasta ->
                reference: tuple(meta_ref, fasta)
                bam_bai: tuple(meta, bam, bai)
            }

        // Run Module
        ARTIC_MAKE_DEPTH_MASK(
            ch_artic_depth_input.bam_bai,
            ch_artic_depth_input.reference
        )
        ch_versions = ch_versions.mix(ARTIC_MAKE_DEPTH_MASK.out.versions.first())
        ch_depth_mask = ARTIC_MAKE_DEPTH_MASK.out.coverage_mask

    } else {
        //
        // PROCESS: Non-Amplicon based data
        //

        // MODULE: Clair3 with no pool splitting
        //
        // Prepare Inputs
        ch_no_pool_input = ch_bam_bai
            .map  { meta, bam, bai -> tuple(meta.ref_id, meta, bam, bai) }
            .combine(ch_ref_with_fai, by: 0)
            .multiMap { _ref_id, meta, bam, bai, meta_ref, fasta , fai ->
                bam_bai: tuple(meta, bam, bai)
                reference: tuple(meta_ref, fasta)
                fai: fai
            }

        // Run Module
        CLAIR3_NO_POOL(
            ch_no_pool_input.bam_bai,
            ch_no_pool_input.reference,
            ch_no_pool_input.fai,
            ch_model
        )
        ch_versions = ch_versions.mix(CLAIR3_NO_POOL.out.versions.first())
        ch_vcf = CLAIR3_NO_POOL.out.vcf


        //
        // MODULE: Make depth mask based on minimum depth to call position
        //
        // Prepare Inputs
        ch_custom_depth_input = ch_bam_bai
            .map  { meta, bam, bai -> tuple(meta.ref_id, meta, bam, bai) }
            .combine(ch_reference.map { meta_ref, fasta -> tuple(meta_ref.id, meta_ref, fasta) }, by: 0)
            .multiMap { _ref_id, meta, bam, bai, meta_ref, fasta ->
                bam_bai: tuple(meta, bam, bai)
                reference: tuple(meta_ref, fasta)
            }

        // Run Module
        CUSTOM_MAKE_DEPTH_MASK(
            ch_custom_depth_input.bam_bai,
            ch_custom_depth_input.reference
        )
        ch_versions = ch_versions.mix(CUSTOM_MAKE_DEPTH_MASK.out.versions.first())
        ch_depth_mask = CUSTOM_MAKE_DEPTH_MASK.out.coverage_mask
        ch_bam = MINIMAP2_ALIGN.out.bam_bai
    }

    //
    // MODULE: Filter VCF variants based on input parameters
    //
    PROCESS_NANOPORE_VCF(
        ch_vcf.join(ch_bam, by:0)
    )
    ch_versions = ch_versions.mix(PROCESS_NANOPORE_VCF.out.versions.first())
    ch_fail_vcf = PROCESS_NANOPORE_VCF.out.fail_vcf

    //
    // MODULE: Mask failing regions
    //
    // Prepare Inputs
    // Join depth and failed VCF
    ch_mask_inputs = ch_depth_mask
        .map { meta, cov -> tuple(meta.id, meta, cov)}
        .join(ch_fail_vcf.map { meta, fail -> tuple(meta.id, fail) })
        .map  { _meta_id, meta, cov, fail -> tuple(meta.ref_id, meta, cov, fail) }

    // Combine with reference
    ch_mask_input = ch_mask_inputs
        .combine(ch_reference.map { meta_ref, fasta -> tuple(meta_ref.id, meta_ref, fasta) }, by: 0)
        .multiMap { _ref_id, meta, cov, fail, meta_ref, fasta ->
            coverage_mask: tuple(meta, cov, fail)
            reference: tuple(meta_ref, fasta)
        }

    // Run Module
    ARTIC_MASK(
        ch_mask_input.coverage_mask,
        ch_mask_input.reference
    )
    ch_versions = ch_versions.mix(ARTIC_MASK.out.versions.first())

    //
    // MODULE: Normalize
    //
    BCFTOOLS_NORM(
        ARTIC_MASK.out.preconsensus
            .join(PROCESS_NANOPORE_VCF.out.pass_vcf, by: [0])
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
    // Prepare Inputs
    ch_adjust_input = BCFTOOLS_CONSENSUS.out.fasta
        .map { meta, con_fasta -> tuple(meta.ref_id, meta, con_fasta) }
        .combine(ch_reference.map { meta_ref, ref_fasta -> tuple(meta_ref.id, meta_ref, ref_fasta) }, by: 0)
        .multiMap { _ref_id, meta, con_fasta, meta_ref, ref_fasta ->
            consensus: tuple(meta, con_fasta)
            reference: tuple(meta_ref, ref_fasta)
        }

    // Run Module
    ADJUST_FASTA_HEADER(
        ch_adjust_input.consensus,
        ch_adjust_input.reference,
        '.consensus',
        ''
    )
    ch_versions = ch_versions.mix(ADJUST_FASTA_HEADER.out.versions.first())

    //
    // MODULE: Create TSV from VCF for final report
    //
    VCF_TO_TSV(
        BCFTOOLS_NORM.out.vcf
            .join(ch_fail_vcf, by: [0])
    )

    emit:
    nanoq_json   = NANOQ.out.stats                      // channel: [ meta, *.json]
    consensus    = ADJUST_FASTA_HEADER.out.consensus    // channel: [ meta, fasta ]
    bam_bai      = ch_bam                               // channel: [ meta, bam, bai ]
    vcf          = BCFTOOLS_NORM.out.vcf                // channel: [ meta, vcf, tbi ]
    variants_tsv = VCF_TO_TSV.out.tsv                   // channel: [ meta, *.tsv ]
    versions     = ch_versions                          // channel: [ path(versions.yml) ]
}
