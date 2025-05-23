//
// Subworkflow for report generation
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO GENERATE REPORTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENERATE_AMPLICON_STATS {

    take:
    ch_bam_bai              //
    ch_consensus            //
    ch_amplicon_bed         //

    main:
    ch_versions = Channel.empty()

    // Setting input up for bedtools coverage
    ch_bam_bai
        .combine(ch_amplicon_bed) // Channel is [meta, bam, bai, bed]
        .map{ it -> [it[0], it[3], it[1]] } // Channel ends [meta, bed, bam]
        .set{ ch_bedtools_input }

    BEDTOOLS_COVERAGE(
        ch_bedtools_input,
        []
    )
    ch_versions = ch_versions.mix(BEDTOOLS_COVERAGE.out.versions.first())

    SAMPLE_AMPLICON_DEPTH(
        BEDTOOLS_COVERAGE.out.bed
    )
    ch_versions = ch_versions.mix(SAMPLE_AMPLICON_DEPTH.out.versions.first())

    AMPLICON_DEPTH_HEATMAP(
        BEDTOOLS_COVERAGE.out.bed.collect{ it[1] }
    )

    SAMPLE_AMPLICON_COMPLETENESS(
        ch_consensus,
        ch_amplicon_bed
    )
    ch_versions = ch_versions.mix(SAMPLE_AMPLICON_COMPLETENESS.out.versions)

    emit:
    versions = ch_versions
}
