/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_measeq_pipeline'
include { SETUP_REFERENCE_DATA    } from '../subworkflows/local/setup_reference_data'
include { NANOPORE_CONSENSUS      } from '../subworkflows/local/nanopore_consensus'
include { ILLUMINA_CONSENSUS      } from '../subworkflows/local/illumina_consensus'
include { FASTQC                  } from '../modules/nf-core/fastqc/main'
include { ADJUST_FASTA_HEADER     } from '../modules/local/artic/subcommands/main'
include { NEXTCLADE_DATASETGET    } from '../modules/nf-core/nextclade/datasetget/main'
include { NEXTCLADE_RUN as NEXTCLADE_RUN_N450   } from '../modules/nf-core/nextclade/run/main'
include { NEXTCLADE_RUN as NEXTCLADE_RUN_CUSTOM } from '../modules/nf-core/nextclade/run/main'
include { SAMTOOLS_DEPTH          } from '../modules/nf-core/samtools/depth/main'
include { COMPARE_INTERNAL_DSID   } from '../modules/local/custom/compare_internal_dsid/main'
include { MAKE_SAMPLE_QC_CSV      } from '../modules/local/qc/sample/main'
include { MAKE_FINAL_QC_CSV       } from '../modules/local/qc/summary/main'
include { GENERATE_AMPLICON_STATS } from '../subworkflows/local/generate_amplicon_stats'
include { GENERATE_REPORT         } from '../subworkflows/local/generate_report'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MEASEQ {

    take:
    ch_input_fastqs // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_metadata = params.metadata ? file(params.metadata, type: 'file', checkIfExists: true) : []
    Channel
        .value(file(params.custom_nextclade_dataset, type: 'dir', checkIfExists: true))
        .set { ch_custom_nextclade_dataset }
    ch_id_fasta = params.dsid_fasta ? file(params.dsid_fasta, type: 'file', checkIfExists: true) : []

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Setup

    //
    // MODULE: Setup nextclade dataset
    //
    NEXTCLADE_DATASETGET(
        params.nextclade_dataset_name,
        params.nextclade_dataset_tag
    )
    ch_versions = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)

    //
    // WORKFLOW: Reference Setup
    //
    SETUP_REFERENCE_DATA(
        params.reference,
        params.primer_bed,
        NEXTCLADE_DATASETGET.out.dataset
    )
    ch_reference                = SETUP_REFERENCE_DATA.out.reference
    ch_reference_fai            = SETUP_REFERENCE_DATA.out.fai
    ch_primer_bed               = SETUP_REFERENCE_DATA.out.primer_bed
    ch_amplicon_bed             = SETUP_REFERENCE_DATA.out.amplicon_bed
    ch_split_amp_pools_bed      = SETUP_REFERENCE_DATA.out.split_amp_pools_bed
    ch_ref_n450                 = SETUP_REFERENCE_DATA.out.ref_n450
    ch_strain                   = SETUP_REFERENCE_DATA.out.strain
    ch_versions                 = ch_versions.mix(SETUP_REFERENCE_DATA.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC(
        ch_input_fastqs
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Consensus Generation

    //
    // WORKFLOWS: Generate consensus and supporting files for either Nanopore or Illumina data
    //
    if( params.platform == 'nanopore' ) {
        //
        // WORKFLOW: Nanopore
        //
        NANOPORE_CONSENSUS(
            ch_reference,
            ch_reference_fai,
            ch_input_fastqs,
            ch_primer_bed,
            ch_split_amp_pools_bed
        )
        ch_read_json   = NANOPORE_CONSENSUS.out.nanoq_json
        ch_bam_bai     = NANOPORE_CONSENSUS.out.bam_bai
        ch_consensus   = NANOPORE_CONSENSUS.out.consensus
        ch_vcf         = NANOPORE_CONSENSUS.out.vcf
        ch_versions    = ch_versions.mix(NANOPORE_CONSENSUS.out.versions)

    } else if( params.platform == 'illumina' ) { 
        //
        // WORKFLOW: Illumina
        //
        ILLUMINA_CONSENSUS(
            ch_reference,
            ch_input_fastqs,
            ch_primer_bed
        )
        ch_read_json   = ILLUMINA_CONSENSUS.out.fastp_json
        ch_bam_bai     = ILLUMINA_CONSENSUS.out.bam_bai
        ch_consensus   = ILLUMINA_CONSENSUS.out.consensus
        ch_vcf         = ILLUMINA_CONSENSUS.out.vcf
        ch_versions    = ch_versions.mix(ILLUMINA_CONSENSUS.out.versions)

    } else {
        error "Please provide the --platform parameter with either 'nanopore' or 'illumina' to run"
    }

    //
    // MODULE: Nextclade Run on generated consensus sequences
    //
    NEXTCLADE_RUN_N450(
        ch_consensus,
        NEXTCLADE_DATASETGET.out.dataset
    )
    ch_versions = ch_versions.mix(NEXTCLADE_RUN_N450.out.versions.first())

    // Using the N450 nextclade align, create renamed, easy to find N450 output
    ADJUST_FASTA_HEADER(
        NEXTCLADE_RUN_N450.out.fasta_aligned,
        ch_reference,
        '.N450',
        '-N450'
    )
    ch_versions = ch_versions.mix(ADJUST_FASTA_HEADER.out.versions.first())

    //
    // MODULE: Run nextclade again using a custom dataset to help determine QC issues in consensus seqs
    //
    NEXTCLADE_RUN_CUSTOM(
        ch_consensus,
        ch_custom_nextclade_dataset
    )
    ch_versions = ch_versions.mix(NEXTCLADE_RUN_CUSTOM.out.versions.first())

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // QC
    SAMTOOLS_DEPTH(
        ch_bam_bai.map{ it -> [it[0], it[1]] },
        [[:], []] // Empty as want to run whole depth
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    //
    // MODULE: Compare to optional internal DSID fasta file to get DSID number
    //
    ch_dsid_tsv = Channel.empty()
    if( params.dsid_fasta ){
        COMPARE_INTERNAL_DSID(
            ADJUST_FASTA_HEADER.out.consensus,
            ch_id_fasta
        )
        ch_dsid_tsv = COMPARE_INTERNAL_DSID.out.tsv
        ch_versions = ch_versions.mix(COMPARE_INTERNAL_DSID.out.versions.first())
    }

    //
    // MODULE: Summarize all of the sample data into 1 CSV file per sample
    //
    MAKE_SAMPLE_QC_CSV(
        ch_bam_bai
            .join(ch_consensus, by: [0])
            .join(SAMTOOLS_DEPTH.out.tsv, by: [0])
            .join(NEXTCLADE_RUN_N450.out.csv, by: [0])
            .join(NEXTCLADE_RUN_CUSTOM.out.csv, by: [0])
            .join(ch_vcf, by: [0])
            .join(ch_read_json, by: [0])
            .join(ch_dsid_tsv, by: [0]).ifEmpty([]),
        ch_strain,
        ch_primer_bed.collect().ifEmpty([])
    )
    ch_versions = ch_versions.mix(MAKE_SAMPLE_QC_CSV.out.versions.first())

    //
    // MODULE: Summarize all individual sample CSVs into 1 final file
    //
    MAKE_FINAL_QC_CSV(
        MAKE_SAMPLE_QC_CSV.out.csv
            .map{ it -> it[1] }
            .collectFile(keepHeader: true, skip: 1, name: 'concat.qc.csv'),
        ch_metadata,
        params.neg_control_pct_threshold,
        params.neg_ctrl_substrings,
        params.skip_negative_grading
    )
    ch_versions = ch_versions.mix(MAKE_FINAL_QC_CSV.out.versions)

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Report Generation

    //
    // WORKFLOW: Amplicon statistics if amplicons were being run
    //
    if( params.primer_bed ){
        GENERATE_AMPLICON_STATS(
            ch_bam_bai,
            ch_consensus,
            ch_amplicon_bed
        )
        ch_versions = ch_versions.mix(GENERATE_AMPLICON_STATS.out.versions)
    }

    //
    // WORKFLOW: Final summary report generation
    //
    GENERATE_REPORT(
        ch_reference,
        ch_reference_fai,
        ch_bam_bai,
        ch_strain,
        SAMTOOLS_DEPTH.out.tsv,
        MAKE_FINAL_QC_CSV.out.csv
    )
    ch_versions = ch_versions.mix(GENERATE_REPORT.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'measeq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
