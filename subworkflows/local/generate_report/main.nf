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
include { softwareVersionsToYAML  } from '../../nf-core/utils_nfcore_pipeline'
include { MAKE_CUSTOM_REPORT      } from '../../../modules/local/custom/custom_report/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO GENERATE REPORTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENERATE_REPORT {

    take:
    ch_reference            // channel: [ meta_ref, fasta ]
    ch_fai                  // channel: [ meta_ref, fai ]
    ch_bam_bai              // channel: [ meta, bam, bai ]
    ch_variants_tsv         // channel: [ meta, tsv ]
    ch_depth_tsv            // channel: [ meta, depth_tsv ]
    ch_overall_qc           // channel: [ csv ]
    ch_versions             // channel: [ version_ymls ]

    main:
    ch_report_template = Channel.fromPath("$projectDir/assets/MeaSeq_Report.Rmd")
    ch_report_subpages = Channel.fromPath("$projectDir/assets/subpage_*.Rmd")

    //
    // MODULE: Calculate the N depth in each position
    //

    // Prepare Inputs
    ch_ref_with_fai = ch_reference
        .map { meta_ref, fasta -> tuple(meta_ref.id, meta_ref, fasta) }
        .join(ch_fai.map { meta_ref, fai -> tuple(meta_ref.id, fai) })
        .map { ref_id, meta_ref, fasta, fai -> tuple(ref_id, meta_ref, fasta, fai) }

    // Combine samples with reference data
    ch_bam_ref = ch_bam_bai
        .map  { meta, bam, bai -> tuple(meta.ref_id, meta, bam, bai) }
        .combine(ch_ref_with_fai, by: 0)
        .multiMap { _ref_id, meta, bam, bai, meta_ref, fasta, fai ->
            ndepth_input: tuple(meta, bam, bai, fasta, fai)
            bam_var_input: tuple(meta, bam, bai, fasta)
        }

    // Run Module
    POSITIONAL_N_DEPTH(
        ch_bam_ref.ndepth_input
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
        ch_bam_ref.bam_var_input
    )
    ch_versions = ch_versions.mix(CALCULATE_BAM_VARIATION.out.versions.first())

    //
    // MODULE: Create a normalized depth matrix using python module
    //
    NORMALIZE_DEPTH_MATRIX(
        ch_depth_tsv.collect{ it[1] }
    )

    //
    // Collate and save software versions
    //  Moved to reporting so that it can be included in the report
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            name:  'measeq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    // Standardize contact information if given
    def website = ''
    if ( params.contact_website ) {
        website = "https://www.${ params.contact_website.replaceFirst(/^https?:\/\//,'').replaceFirst(/^www\./,'') }"
    }
    def email = ''
    if (params.contact_email) {
        email = params.contact_email ==~ /[^@]+@[^@]+\.[^@]+/ ? "mailto:${params.contact_email}" : ""
    }

    //
    // MODULE: Make custom final HTML Report
    //
    revision = workflow.revision ? workflow.revision : 'main'
    MAKE_CUSTOM_REPORT(
        ch_overall_qc,
        ch_depth_tsv.collect{ it[1] },
        POSITIONAL_N_DEPTH.out.tsv.collect{ it[1] },
        PYSAMSTATS.out.tsv.collect{ it[1] },
        CALCULATE_BAM_VARIATION.out.csv.collect{ it[1] },
        ch_variants_tsv.collect{ it[1] },
        NORMALIZE_DEPTH_MATRIX.out.full_csv,
        ch_report_template,
        ch_report_subpages.collect(),
        ch_collated_versions,
        workflow.manifest.version,
        revision,
        nextflow.version,
        params.contact_name ?: '',
        email,
        params.contact_phone ?: '',
        website
    )

    emit:
    versions = ch_versions      // channel: [ path(versions.yml) ]
}
