//
// Subworkflow for report generation
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PYSAMSTATS         } from '../../../modules/local/pysamstats/main'
include { POSITIONAL_N_DEPTH } from '../../../modules/local/custom/positional_n_depth/main'
include { MAKE_CUSTOM_REPORT } from '../../../modules/local/custom/custom_report/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO GENERATE REPORTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENERATE_REPORT {

    take:
    ch_reference             //
    ch_reference_fai         //
    ch_bam_bai               //
    ch_strain                //
    ch_depth_tsv             //
    ch_overall_qc            //

    main:
    ch_versions = Channel.empty()
    ch_report_template = Channel.fromPath("$projectDir/assets/MeaSeq_Report.Rmd")
    ch_report_subpages = Channel.fromPath("$projectDir/assets/subpage_*.Rmd")

    POSITIONAL_N_DEPTH(
        ch_bam_bai,
        ch_reference,
        ch_reference_fai
    )
    ch_versions = ch_versions.mix(POSITIONAL_N_DEPTH.out.versions.first())

    PYSAMSTATS(
        ch_bam_bai,
        'baseq'
    )
    ch_versions = ch_versions.mix(PYSAMSTATS.out.versions.first())

    // Report is broken for the moment which is expected with the changes
    // MAKE_CUSTOM_REPORT(
    //     ch_overall_qc,
    //     ch_depth_tsv.collect{ it[1] },
    //     POSITIONAL_N_DEPTH.out.tsv.collect{ it[1] },
    //     PYSAMSTATS.out.tsv.collect{ it[1] },
    //     ch_strain,
    //     ch_report_template,
    //     ch_report_subpages.collect()
    // )
    // ch_versions = ch_versions.mix(MAKE_CUSTOM_REPORT.out.versions)

    emit:
    versions = ch_versions
}
