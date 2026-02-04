//
// Subworkflow for amplicon report generation
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BEDTOOLS_COVERAGE             } from '../../../modules/nf-core/bedtools/coverage/main'
include { SAMPLE_AMPLICON_DEPTH         } from '../../../modules/local/custom/sample_amplicon_depth/main'
include { AMPLICON_DEPTH_HEATMAP        } from '../../../modules/local/custom/amplicon_depth_heatmap/main'
include { SAMPLE_AMPLICON_COMPLETENESS  } from '../../../modules/local/custom/sample_amplicon_completeness/main'
include { AMPLICON_COMPLETENESS_HEATMAP } from '../../../modules/local/custom/amplicon_completeness_heatmap/main'
include { AMPLICON_MULTIQC              } from '../../../modules/local/multiqc/amplicon_multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO GENERATE REPORTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENERATE_AMPLICON_STATS {

    take:
    ch_bam_bai              // channel: [ meta, bam, bai ]
    ch_consensus            // channel: [ meta, fasta ]
    ch_amplicon_bed         // channel: [ meta_ref, bed ]

    main:
    ch_versions         = Channel.empty()
    ch_multiqc_files    = Channel.empty()
    ch_multiqc_config   = Channel.fromPath("$projectDir/assets/amplicon_multiqc_config.yml", checkIfExists: true)

    //
    // MODULE: Run bedtools coverage
    //
    // Prepare Inputs
    ch_bedtools_input = ch_bam_bai
        .map{ meta, bam, _bai -> tuple(meta.ref_id, meta, bam) }
        .combine(ch_amplicon_bed.map { meta_ref, bed -> tuple(meta_ref.id, meta_ref, bed) }, by: 0)
        .map { _ref_id, meta, bam, _meta_ref, bed ->
            tuple(meta, bed, bam)
        }

    // Run Module
    BEDTOOLS_COVERAGE(
        ch_bedtools_input,
        []
    )
    ch_versions = ch_versions.mix(BEDTOOLS_COVERAGE.out.versions.first())

    //
    // MODULE: Calculate the per amplicon depth for each sample with csvtk
    //
    SAMPLE_AMPLICON_DEPTH(
        BEDTOOLS_COVERAGE.out.bed
    )
    ch_versions = ch_versions.mix(SAMPLE_AMPLICON_DEPTH.out.versions.first())
    ch_depth_tsvs = SAMPLE_AMPLICON_DEPTH.out.tsv
        .map { meta, tsv -> tuple(meta.ref_id, tsv) }
        .groupTuple()

    //
    // MODULE: Summarize and reformat amplicon depth into a matrix for multiqc heatmap reporting
    //
    // Prepare Inputs
    ch_depth_heatmap_input = BEDTOOLS_COVERAGE.out.bed
        .map { meta, path -> tuple(meta.ref_id, path) }
        .groupTuple()

    // Run Module
    AMPLICON_DEPTH_HEATMAP(
        ch_depth_heatmap_input
    )

    //
    // MODULE: Calculate the per-amplicon completeness with custom python script
    //
    // Prepare Inputs
    ch_completeness_input = ch_consensus
        .map { meta, con_fasta -> tuple(meta.ref_id, meta, con_fasta) }
        .combine(ch_amplicon_bed.map { meta_ref, bed -> tuple(meta_ref.id, bed) }, by: 0)
        .map { _ref_id, meta, con_fasta, bed -> tuple(meta, con_fasta, bed) }

    // Run Module
    SAMPLE_AMPLICON_COMPLETENESS(
        ch_completeness_input
    )
    ch_versions = ch_versions.mix(SAMPLE_AMPLICON_COMPLETENESS.out.versions)
    ch_completeness_tsvs = SAMPLE_AMPLICON_COMPLETENESS.out.tsv
        .map { meta, tsv -> tuple(meta.ref_id, tsv) }
        .groupTuple()

    //
    // MODULE: Summarize the per-amplicon completeness with csvtk into a matrix for multiqc heatmap
    //
    // Prepare Inputs
    ch_completeness_heatmap_input = SAMPLE_AMPLICON_COMPLETENESS.out.tsv
        .map { meta, tsv -> tuple(meta.ref_id, tsv) }
        .groupTuple()

    // Run Module
    AMPLICON_COMPLETENESS_HEATMAP(
        ch_completeness_heatmap_input
    )

    //
    // MODULE: Run amplicon multiqc for reporting
    //
    // Prepare Inputs
    ch_amplicon_multiqc_input = ch_depth_tsvs
        .join(AMPLICON_DEPTH_HEATMAP.out.heatmap_tsv)
        .join(ch_completeness_tsvs)
        .join(AMPLICON_COMPLETENESS_HEATMAP.out.heatmap_tsv)
        .map{ ref_id, depth_tsvs, depth_heatmap, completeness_tsvs, completeness_heatmap ->
            tuple(ref_id, [depth_tsvs, depth_heatmap, completeness_tsvs, completeness_heatmap].flatten())
        }

    // Run Module
    AMPLICON_MULTIQC(
        ch_amplicon_multiqc_input,
        ch_multiqc_config.collect()
    )

    emit:
    versions = ch_versions  // channel: [ path(versions.yml) ]
}
