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
    ch_bam_bai              // channel: [ [id], bam, bai ]
    ch_consensus            // channel: [ [id], fasta ]
    ch_amplicon_bed         // channel: [ bed ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/amplicon_multiqc_config.yml", checkIfExists: true)

    // Setting input up for bedtools coverage
    ch_bam_bai
        .combine(ch_amplicon_bed) // Channel is [meta, bam, bai, bed]
        .map{ it -> [it[0], it[3], it[1]] } // Channel ends [meta, bed, bam]
        .set{ ch_bedtools_input }

    //
    // MODULE: Run bedtools coverage
    //
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
    ch_multiqc_files = ch_multiqc_files.mix(SAMPLE_AMPLICON_DEPTH.out.tsv.collect{ it -> it[1] })

    //
    // MODULE: Summarize and reformat amplicon depth into a matrix for multiqc heatmap reporting
    //
    AMPLICON_DEPTH_HEATMAP(
        BEDTOOLS_COVERAGE.out.bed.collect{ it[1] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(AMPLICON_DEPTH_HEATMAP.out.heatmap_tsv)

    //
    // MODULE: Calculate the per-amplicon completeness with custom python script
    //
    SAMPLE_AMPLICON_COMPLETENESS(
        ch_consensus,
        ch_amplicon_bed
    )
    ch_versions = ch_versions.mix(SAMPLE_AMPLICON_COMPLETENESS.out.versions)

    //
    // MODULE: Summarize the per-amplicon completeness with csvtk into a matrix for multiqc heatmap
    //
    AMPLICON_COMPLETENESS_HEATMAP(
        SAMPLE_AMPLICON_COMPLETENESS.out.tsv.collect{ it -> it[1] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(AMPLICON_COMPLETENESS_HEATMAP.out.heatmap_tsv)

    //
    // MODULE: Run amplicon multiqc for reporting
    //
    AMPLICON_MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config
    )

    emit:
    versions = ch_versions
}
