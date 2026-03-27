/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SETUP_REFERENCE_DATA    } from '../subworkflows/local/setup_reference_data'
include { NANOPORE_CONSENSUS      } from '../subworkflows/local/nanopore_consensus'
include { ILLUMINA_CONSENSUS      } from '../subworkflows/local/illumina_consensus'
include { FASTQC                  } from '../modules/nf-core/fastqc/main'
include { ADJUST_FASTA_HEADER as ADJUST_N450_FASTA_HEADER } from '../modules/local/artic/subcommands/main'
include { NEXTCLADE_DATASETGET    } from '../modules/nf-core/nextclade/datasetget/main'
include { NEXTCLADE_RUN as NEXTCLADE_RUN_N450             } from '../modules/nf-core/nextclade/run/main'
include { NEXTCLADE_RUN as NEXTCLADE_RUN_CUSTOM           } from '../modules/nf-core/nextclade/run/main'
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
    samplesheet   // channel: samplesheet read in from --input // [ meta(id, single-end), fastqs(f1,f2) ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_metadata = params.metadata ? file(params.metadata, type: 'file', checkIfExists: true) : []
    ch_custom_nextclade_dataset = Channel
        .value(file("$projectDir/assets/custom_measles_nextclade_dataset", type: 'dir', checkIfExists: true))
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
        samplesheet,
        NEXTCLADE_DATASETGET.out.dataset
    )
    ch_samples              = SETUP_REFERENCE_DATA.out.samples
    ch_primer_bed           = SETUP_REFERENCE_DATA.out.primer_bed
    ch_reference            = SETUP_REFERENCE_DATA.out.reference
    ch_fai                  = SETUP_REFERENCE_DATA.out.fai
    ch_refstats             = SETUP_REFERENCE_DATA.out.refstats
    ch_genome_bed           = SETUP_REFERENCE_DATA.out.genome_bed
    ch_amplicon_bed         = SETUP_REFERENCE_DATA.out.amplicon_bed
    ch_split_amp_pools_bed  = SETUP_REFERENCE_DATA.out.split_amp_pools_bed
    ch_ref_n450             = SETUP_REFERENCE_DATA.out.ref_n450
    ch_versions             = ch_versions.mix(SETUP_REFERENCE_DATA.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC(
        ch_samples
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
            ch_fai,
            ch_samples,
            ch_primer_bed,
            ch_split_amp_pools_bed
        )
        ch_read_json    = NANOPORE_CONSENSUS.out.nanoq_json
        ch_bam_bai      = NANOPORE_CONSENSUS.out.bam_bai
        ch_consensus    = NANOPORE_CONSENSUS.out.consensus
        ch_vcf          = NANOPORE_CONSENSUS.out.vcf
        ch_variants_tsv = NANOPORE_CONSENSUS.out.variants_tsv
        ch_versions     = ch_versions.mix(NANOPORE_CONSENSUS.out.versions)

    } else if( params.platform == 'illumina' ) {
        //
        // WORKFLOW: Illumina
        //
        ILLUMINA_CONSENSUS(
            ch_reference,
            ch_fai,
            ch_samples,
            ch_primer_bed
        )
        ch_read_json    = ILLUMINA_CONSENSUS.out.fastp_json
        ch_bam_bai      = ILLUMINA_CONSENSUS.out.bam_bai
        ch_consensus    = ILLUMINA_CONSENSUS.out.consensus
        ch_vcf          = ILLUMINA_CONSENSUS.out.vcf
        ch_variants_tsv = ILLUMINA_CONSENSUS.out.variants_tsv
        ch_versions     = ch_versions.mix(ILLUMINA_CONSENSUS.out.versions)

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

    //
    // MODULE: Adjust N450 sequence headers to make downstream processes easier
    //

    // Prepare Input
    // Using the N450 nextclade align, create renamed & easy to find N450 output
    ch_n450_adjust_input = NEXTCLADE_RUN_N450.out.fasta_aligned
        .map { meta, n450_fasta -> tuple(meta.ref_id, meta, n450_fasta) }
        .combine(ch_reference.map { meta_ref, ref_fasta -> tuple(meta_ref.id, meta_ref, ref_fasta) }, by: 0)
        .multiMap { _ref_id, meta, n450_fasta, meta_ref, ref_fasta ->
            n450: tuple(meta, n450_fasta)
            reference: tuple(meta_ref, ref_fasta)
        }

    // Run Module
    ADJUST_N450_FASTA_HEADER(
        ch_n450_adjust_input.n450,
        ch_n450_adjust_input.reference,
        '.N450',
        '-N450'
    )
    ch_versions = ch_versions.mix(ADJUST_N450_FASTA_HEADER.out.versions.first())

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

    //
    // MODULE: get the sequencing depth of each sample
    SAMTOOLS_DEPTH(
        ch_bam_bai.map{ meta, bam, _bai -> tuple(meta, bam) },
        [[:], []] // Empty as we want to run whole depth
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    //
    // MODULE: Compare to optional internal DSID fasta file to get DSID number
    //
    COMPARE_INTERNAL_DSID(
        ADJUST_N450_FASTA_HEADER.out.consensus
            .map{ it -> it[1] }
            .collectFile(name: 'N450.fasta', sort: { it.baseName }),
        ch_id_fasta
    )
    ch_dsid_results = COMPARE_INTERNAL_DSID.out.dsid_tsv
    ch_versions = ch_versions.mix(COMPARE_INTERNAL_DSID.out.versions.first())

    //
    // MODULE: Summarize all of the sample data into 1 CSV file per sample
    //

    //Prepare Inputs
    ch_overall_sample = ch_bam_bai
        .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai)}
        .join(ch_consensus.map { meta, con_fasta -> tuple(meta.id, con_fasta)}, by: [0])
        .join(ADJUST_N450_FASTA_HEADER.out.consensus.map { meta, n450_fasta -> tuple(meta.id, n450_fasta)}, by: [0])
        .join(SAMTOOLS_DEPTH.out.tsv.map { meta, depth_bed -> tuple(meta.id, depth_bed)}, by: [0])
        .join(NEXTCLADE_RUN_N450.out.csv.map { meta, nextclade_n450 -> tuple(meta.id, nextclade_n450)}, by: [0])
        .join(NEXTCLADE_RUN_CUSTOM.out.csv.map { meta, nextclade_full -> tuple(meta.id, nextclade_full)}, by: [0])
        .join(ch_vcf.map { meta, vcf, tbi -> tuple(meta.id, vcf, tbi)}, by: [0])
        .join(ch_read_json.map { meta, read_json -> tuple(meta.id, read_json)}, by: [0])
        .join(ch_variants_tsv.map { meta, tsv -> tuple(meta.id, tsv)}, by: [0])
        .map { _meta_id, meta, bam, bai, con_fasta, n450_fasta, depth_bed, nextclade_n450, nextclade_full, vcf, tbi, read_json, var_tsv ->
            tuple(meta.ref_id, meta, bam, bai, con_fasta, n450_fasta, depth_bed, nextclade_n450, nextclade_full, vcf, tbi, read_json, var_tsv)
        }

    // Prepare primers file and join with genotype if available or just output genotype
    if ( params.amplicon || params.primer_bed ) {
        ch_overall_ref = ch_reference
            .map{ meta_ref, _fasta -> tuple(meta_ref.id, meta_ref.genotype) }
            .join(ch_primer_bed.map { meta_ref, bed -> tuple(meta_ref.id, bed)}, by: [0])
            .map { ref_id, ref_genotype, bed ->
                tuple(ref_id, ref_genotype, bed)
            }
            .unique { it[0] }
    } else {
        ch_overall_ref = ch_reference
            .map{ meta_ref, _fasta -> tuple(meta_ref.id, meta_ref.genotype, []) }
    }

    // Combine genotype/primers with their corresponding samples according to ref_id
    ch_sample_qc_input = ch_overall_sample
        .combine(ch_overall_ref, by: 0)
        .map { _ref_id, meta, bam, bai, con_fasta, n450_fasta, depth_bed, nextclade_n450, nextclade_full, vcf, tbi, read_json, var_tsv, genotype, bed ->
            tuple(meta, bam, bai, con_fasta, n450_fasta, depth_bed, nextclade_n450, nextclade_full, vcf, tbi, read_json, var_tsv, genotype, bed)
        }

    // Run Module
    MAKE_SAMPLE_QC_CSV(
        ch_sample_qc_input,
        ch_dsid_results.collect()
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
    if( params.amplicon || params.primer_bed ) {
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
        ch_fai,
        ch_bam_bai,
        ch_variants_tsv,
        SAMTOOLS_DEPTH.out.tsv,
        MAKE_FINAL_QC_CSV.out.csv,
        ch_versions
    )
    ch_versions = ch_versions.mix(GENERATE_REPORT.out.versions)

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
