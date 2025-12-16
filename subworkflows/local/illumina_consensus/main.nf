//
// Subworkflow for amplicon and non-amplicon consensus sequence generation for Illumina data
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Initial Steps
include { FASTP                     } from '../../../modules/nf-core/fastp/main'
// Amplicon Specific
include { IVAR_TRIM                 } from '../../../modules/local/ivar/trim/main'
// Variant Calling and Consensus Generation
include { BWAMEM2_INDEX             } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM               } from '../../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_SORT             } from '../../../modules/nf-core/samtools/sort/main'
include { BAM_MARKDUPLICATES_PICARD } from '../../../subworkflows/nf-core/bam_markduplicates_picard/main'
include { BAM_STATS_SAMTOOLS        } from '../../../subworkflows/nf-core/bam_stats_samtools/main'
include { FREEBAYES                 } from '../../../modules/local/freebayes/main'
include { PROCESS_VCF               } from '../../../modules/local/process_vcf/main'
include { CUSTOM_MAKE_DEPTH_MASK    } from '../../../modules/local/artic/subcommands/main'
include { BCFTOOLS_CONSENSUS as BCFTOOLS_CONSENSUS_AMBIGUOUS } from '../../../modules/nf-core/bcftools/consensus/main'
include { BCFTOOLS_CONSENSUS as BCFTOOLS_CONSENSUS_FINAL     } from '../../../modules/nf-core/bcftools/consensus/main'
include { ADJUST_FASTA_HEADER       } from '../../../modules/local/artic/subcommands/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO SETUP REFERENCE DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ILLUMINA_CONSENSUS {

    take:
    ch_reference    // channel: [ meta_ref, fasta ]
    ch_fai          // channel: [ meta_ref, fai ]
    ch_samples      // channel: [ meta[id, single_end, ref_id], fastqs ]
    ch_primer_bed   // channel: [ meta_ref, bed ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run fastp for read quality filtering
    //
    FASTP(
        ch_samples,
        [],
        '',
        '',
        ''
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    //
    // MODULE: Create BWAMEM index of reference
    //
    BWAMEM2_INDEX(
        ch_reference
    )
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

    //
    // MODULE: Run BWAMEM to map to reference
    //  Using modules.config to filter reads with view
    //
    // Prepare Input
    // Link FASTP outputs with INDEX & Join reads to their matching reference by ref_id
    ch_reads_by_ref = FASTP.out.reads
        .map { meta, reads -> tuple(meta.ref_id, meta, reads) }

    ch_ref_with_index = ch_reference
        .map { meta_ref, fasta -> tuple(meta_ref.id, meta_ref, fasta) }
        .join(
            BWAMEM2_INDEX.out.index.map { meta_ref, index_dir -> tuple(meta_ref.id, index_dir) }
        )
        .map { ref_id, meta_ref, fasta, index_dir ->
            tuple(ref_id, meta_ref, fasta, index_dir)
        }

    // Create BWAMEM2 expected input
    ch_bwa_mem_inputs = ch_reads_by_ref
        .combine(ch_ref_with_index, by: 0)
        .multiMap { ref_id, meta, reads, meta_ref, fasta, index_dir ->
            reads: tuple(meta, reads)
            index: tuple(meta_ref, index_dir)
            reference: tuple(meta_ref, fasta)
            mode : ''
        }


    // Run Module
    BWAMEM2_MEM(
        ch_bwa_mem_inputs.reads,
        ch_bwa_mem_inputs.index,
        ch_bwa_mem_inputs.reference,
        ch_bwa_mem_inputs.mode
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

    //
    // MODULE: Sort and index output bam file from BWA
    //
    // Prepare Inputs
    ch_samtools_sort_input = BWAMEM2_MEM.out.bam
        .map { meta, bam -> tuple(meta.ref_id, meta, bam) }
        .combine(
            ch_reference.map { meta_ref, fasta -> tuple(meta_ref.id, meta_ref, fasta) }, by: [0]
        )
        .multiMap { _ref_id, meta, bam, meta_ref, fasta ->
            bam: tuple(meta, bam)
            reference: tuple(meta_ref, fasta)
            bai: 'bai'
        }

    // Run Module
    SAMTOOLS_SORT(
        ch_samtools_sort_input.bam,
        ch_samtools_sort_input.reference,
        ch_samtools_sort_input.bai
    )
    ch_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_SORT.out.bai, by: [0])
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    //
    // PROCESS: if we have amplicon data, want to make sure that the primers are not affecting
    //  the output results so trimming them is required
    //
    if ( params.amplicon || params.primer_bed ) {
        //
        // MODULE: Run IVAR Trim
        //
        // Prepare Inputs
        ch_ivar_trim_input = ch_bam_bai
            .map { meta, bam, bai -> tuple(meta.ref_id, meta, bam, bai) }
            .combine(ch_primer_bed.map { meta_ref, bed -> tuple(meta_ref.id, bed) }, by: 0)
            .map { _ref_id, meta, bam, bai, bed ->
                tuple(meta, bam, bai, bed)
            }

        // Run Module
        IVAR_TRIM(
            ch_ivar_trim_input
        )
        ch_bam_bai = IVAR_TRIM.out.bam
        ch_versions = ch_versions.mix(IVAR_TRIM.out.versions)
    }

    //
    // SUBWORKFLOW: Mark duplicates if arg is given
    //
    // Prepare Input
    ch_bam_bai_by_ref = ch_bam_bai
        .map { meta, bam, bai -> tuple(meta.ref_id, meta, bam, bai) }

    // Join reference index file with fasta file
    ch_ref_with_fai = ch_reference
        .map { meta_ref, fasta -> tuple(meta_ref.id, meta_ref, fasta) }
        .join(ch_fai.map { meta_ref, fai -> tuple(meta_ref.id, fai) }, by: [0])
        .map { ref_id, meta_ref, fasta, fai ->
                tuple(ref_id, meta_ref, fasta, fai)
        }

    // Create expected inputs for the modules
    ch_bam_input = ch_bam_bai_by_ref
        .combine(ch_ref_with_fai, by: 0)
        .multiMap { _ref_id, meta, bam, bai, meta_ref, fasta, fai ->
            bam  : tuple(meta, bam)
            reference: tuple(meta_ref, fasta)
            fai: tuple(meta_ref, fai)
            bam_bai: tuple(meta, bam, bai)
        }

    if( params.remove_duplicates ) {
        // Run Subworkflow
        BAM_MARKDUPLICATES_PICARD(
            ch_bam_input.bam,
            ch_bam_input.reference,
            ch_bam_input.fai
        )
        ch_bam_bai = BAM_MARKDUPLICATES_PICARD.out.bam
                        .join(BAM_MARKDUPLICATES_PICARD.out.bai, by: [0])
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)
    } else {
        // Run subworkflow: Get samtools stats
        //  These are also in the PICARD workflow
        BAM_STATS_SAMTOOLS(
            ch_bam_input.bam_bai,
            ch_bam_input.reference
        )
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)
    }

    //
    // MODULE: Run Freebayes to call variants
    //
    // Prepare Input for both FREEBAYES & CUSTOM_MAKE_DEPTH_MASK modules
    ch_bam_bai_by_ref = ch_bam_bai
        .map { meta, bam, bai -> tuple(meta.ref_id, meta, bam, bai) }
    ch_freebayes_depth_input = ch_bam_bai_by_ref
        .combine(ch_ref_with_fai, by:0)
        .multiMap { _ref_id, meta, bam, bai, meta_ref, fasta, fai ->
            reference: tuple(meta_ref, fasta)
            bam_bai: tuple(meta, bam, bai)
        }

    // Run Module
    FREEBAYES(
        ch_freebayes_depth_input.bam_bai,
        ch_freebayes_depth_input.reference
    )
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)

    //
    // MODULE: Process freebayes variant calls with custom python script and bcftools norm
    //
    // Prepare Input
    ch_process_vcf_input = FREEBAYES.out.vcf
            .map { meta, vcf -> tuple(meta.ref_id, meta, vcf) }
            .combine(ch_reference.map { meta_ref, fasta -> tuple(meta_ref.id, meta_ref, fasta) }, by: 0)
            .multiMap { _ref_id, meta, vcf, meta_ref, fasta ->
                vcf: tuple(meta, vcf)
                reference: tuple(meta_ref, fasta)
            }

    // Run Module
    PROCESS_VCF(
        ch_process_vcf_input.vcf,
        ch_process_vcf_input.reference
    )
    ch_versions = ch_versions.mix(PROCESS_VCF.out.versions)

    //
    // MODULE: Make a depth mask based on the minimum depth required to call a position
    //
    CUSTOM_MAKE_DEPTH_MASK(
        ch_freebayes_depth_input.bam_bai,
        ch_freebayes_depth_input.reference
    )
    ch_versions = ch_versions.mix(CUSTOM_MAKE_DEPTH_MASK.out.versions)

    //
    // MODULE: Create intermediate fasta file with IUPACs for ambiguous positions from freebayes
    //
    // Prepare Input
    ch_ambiguous_vcf_restructured = PROCESS_VCF.out.ambiguous_vcf
            .map {meta, vcf_gz, tbi -> tuple(meta.ref_id, meta, vcf_gz, tbi) }
            .combine(ch_reference.map { meta_ref, fasta -> tuple(meta_ref.id, meta_ref, fasta) }, by: 0)
            .map { _ref_id, meta, vcf_gz, tbi, _meta_ref, fasta ->
                tuple(meta, vcf_gz, tbi, fasta, [])
            }

    // Run Module
    BCFTOOLS_CONSENSUS_AMBIGUOUS(
        ch_ambiguous_vcf_restructured
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS_AMBIGUOUS.out.versions)

    //
    // MODULE: Create final consensus sequence with all variants
    //
    BCFTOOLS_CONSENSUS_FINAL(
        PROCESS_VCF.out.consensus_vcf
            .join(BCFTOOLS_CONSENSUS_AMBIGUOUS.out.fasta, by: [0])
            .join(CUSTOM_MAKE_DEPTH_MASK.out.coverage_mask, by: [0])
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS_FINAL.out.versions)

    //
    // MODULE: Adjust final consensus sequence headers to make downstream processes easier
    //
    // Prepare Input
    ch_adjust_input = BCFTOOLS_CONSENSUS_FINAL.out.fasta
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
    ch_versions = ch_versions.mix(ADJUST_FASTA_HEADER.out.versions)

    emit:
    fastp_json   = FASTP.out.json                       // channel: [ meta, *.json]
    bam_bai      = ch_bam_bai                           // channel: [ meta, bam, bai ]
    consensus    = ADJUST_FASTA_HEADER.out.consensus    // channel: [ meta, fasta ]
    vcf          = PROCESS_VCF.out.total_vcf            // channel: [ meta, vcf, tbi ]
    variants_tsv = PROCESS_VCF.out.variants_tsv         // channel: [ meta, *.tsv ]
    versions     = ch_versions                          // channel: [ path(versions.yml) ]
}
