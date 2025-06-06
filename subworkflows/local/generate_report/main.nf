//
// Subworkflow for custom report generation
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PYSAMSTATS              } from '../../../modules/local/pysamstats/main'
include { POSITIONAL_N_DEPTH      } from '../../../modules/local/custom/positional_n_depth/main'
include { CALCULATE_BAM_VARIATION } from '../../../modules/local/custom/calculate_bam_variation/main'
include { NORMALIZE_DEPTH_MATRIX  } from '../../../modules/local/custom/normalize_depth_matrix/main'
include { MAKE_CUSTOM_REPORT      } from '../../../modules/local/custom/custom_report/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO GENERATE REPORTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENERATE_REPORT {

    take:
    ch_reference            // channel: [ [id], fasta ]
    ch_reference_fai        // channel: [ fai ]
    ch_bam_bai              // channel: [ [id], bam, bai ]
    ch_variants_tsv         // channel: [ [id], tsv ]
    ch_genotype             // channel: [ genotype ]
    ch_depth_tsv            // channel: [ [id], depth_tsv ]
    ch_overall_qc           // channel: [ csv ]

    main:
    ch_versions = Channel.empty()
    ch_report_template = Channel.fromPath("$projectDir/assets/MeaSeq_Report.Rmd")
    ch_report_subpages = Channel.fromPath("$projectDir/assets/subpage_*.Rmd")

    //
    // MODULE: Calculate the N depth in each position
    //
    POSITIONAL_N_DEPTH(
        ch_bam_bai,
        ch_reference,
        ch_reference_fai
    )
    ch_versions = ch_versions.mix(POSITIONAL_N_DEPTH.out.versions.first())

    //
    // MODULE: Calculate the per-base quality from the bam file
    //
    PYSAMSTATS(
        ch_bam_bai,
        'baseq'
    )
    ch_versions = ch_versions.mix(PYSAMSTATS.out.versions.first())

    //
    // MODULE: Calculate the underlying variation in the BAM file
    //
    CALCULATE_BAM_VARIATION(
        ch_bam_bai,
        ch_reference
    )
    ch_versions = ch_versions.mix(CALCULATE_BAM_VARIATION.out.versions.first())

    //
    // MODULE: Create a normalized depth matrix using python module
    //
    NORMALIZE_DEPTH_MATRIX(
        ch_depth_tsv.collect{ it[1] }
    )

    //
    // MODULE: Make custom final HTML Report
    //
    MAKE_CUSTOM_REPORT(
        ch_overall_qc,
        ch_depth_tsv.collect{ it[1] },
        POSITIONAL_N_DEPTH.out.tsv.collect{ it[1] },
        PYSAMSTATS.out.tsv.collect{ it[1] },
        CALCULATE_BAM_VARIATION.out.csv.collect{ it[1] },
        ch_variants_tsv.collect{ it[1] },
        NORMALIZE_DEPTH_MATRIX.out.full_csv,
        ch_genotype,
        ch_report_template,
        ch_report_subpages.collect()
    )
    ch_versions = ch_versions.mix(MAKE_CUSTOM_REPORT.out.versions)

    emit:
    versions = ch_versions
}
