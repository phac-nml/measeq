//
// Subworkflow for amplicon and non-amplicon consensus sequence generation for Illumina data
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Initial Steps
include { FASTP                   } from '../../../modules/nf-core/fastp/main'
// Amplicon Specific
include { IVAR_TRIM               } from '../../../modules/local/ivar/trim/main'
// Variant Calling and Consensus Generation
include { BWAMEM2_INDEX           } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM             } from '../../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_INDEX          } from '../../../modules/nf-core/samtools/index/main'
include { FREEBAYES               } from '../../../modules/local/freebayes/main'
include { PROCESS_VCF             } from '../../../modules/local/process_vcf/main'
include { CUSTOM_MAKE_DEPTH_MASK  } from '../../../modules/local/artic/subcommands/main'
include { BCFTOOLS_CONSENSUS as BCFTOOLS_CONSENSUS_AMBIGUOUS } from '../../../modules/nf-core/bcftools/consensus/main'
include { BCFTOOLS_CONSENSUS as BCFTOOLS_CONSENSUS_FINAL     } from '../../../modules/nf-core/bcftools/consensus/main'
include { ADJUST_FASTA_HEADER     } from '../../../modules/local/artic/subcommands/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO SETUP REFERENCE DATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ILLUMINA_CONSENSUS {

    take:
    ch_reference           // channel: [ [id], fasta ]
    ch_input_fastqs        // channel: [ [id, single_end], fastqs ]
    ch_primer_bed          // channel: [ bed ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run fastp for read quality filtering
    //
    FASTP(
        ch_input_fastqs,
        [],
        '',
        '',
        ''
    )

    //
    // MODULE: Create BWAMEM index of reference
    //
    BWAMEM2_INDEX(
        ch_reference
    )

    //
    // MODULE: Run BWAMEM to map to reference
    //
    BWAMEM2_MEM(
        FASTP.out.reads,
        BWAMEM2_INDEX.out.index,
        ch_reference,
        'sort'
    )

    //
    // MODULE: Index output bam file from BWA
    //
    SAMTOOLS_INDEX(
        BWAMEM2_MEM.out.bam
    )
    ch_bam_bai = BWAMEM2_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0])

    //
    // PROCESS: if we have amplicon data, want to make sure that the primers are not affecting
    //  the output results so trimming them is required
    //
    if( params.primer_bed ) {
        // iVar Trim
        IVAR_TRIM(
            ch_bam_bai,
            ch_primer_bed
        )
        ch_bam_bai = IVAR_TRIM.out.bam
    }

    //
    // MODULE: Run Freebayes to call variants
    //
    FREEBAYES(
        ch_bam_bai,
        ch_reference
    )

    //
    // MODULE: Process freebayes variant calls with custom python script and bcftools norm
    //
    PROCESS_VCF(
        FREEBAYES.out.vcf,
        ch_reference
    )

    // For some reason it did not want to let me join the [] alone to the new channel for the next steps
    //  So we've done it this way and that works
    PROCESS_VCF.out.ambiguous_vcf
        .combine(ch_reference.map{it -> it[1]})
        .map{ it -> [it[0], it[1], it[2], it[3], []]}
        .set{ ch_ambiguous_vcf_restructured }

    //
    // MODULE: Make a depth mask based on the minimum depth required to call a position
    //
    CUSTOM_MAKE_DEPTH_MASK(
        ch_bam_bai,
        ch_reference
    )

    //
    // MODULE: Create intermediate fasta file with IUPACs for ambiguous positions from freebayes
    //
    BCFTOOLS_CONSENSUS_AMBIGUOUS(
        ch_ambiguous_vcf_restructured
    )

    //
    // MODULE: Create final consensus sequence with all variants
    //
    BCFTOOLS_CONSENSUS_FINAL(
        PROCESS_VCF.out.consensus_vcf
            .join(BCFTOOLS_CONSENSUS_AMBIGUOUS.out.fasta, by: [0])
            .join(CUSTOM_MAKE_DEPTH_MASK.out.coverage_mask, by: [0])
    )

    //
    // MODULE: Adjust final consensus sequence headers to make downstream processes easier
    //
    ADJUST_FASTA_HEADER(
        BCFTOOLS_CONSENSUS_FINAL.out.fasta,
        ch_reference,
        '.consensus',
        ''
    )

    emit:
    fastp_json   = FASTP.out.json
    bam_bai      = ch_bam_bai
    consensus    = ADJUST_FASTA_HEADER.out.consensus
    vcf          = PROCESS_VCF.out.consensus_vcf
    variants_tsv = PROCESS_VCF.out.variants_tsv
    versions     = ch_versions
}
