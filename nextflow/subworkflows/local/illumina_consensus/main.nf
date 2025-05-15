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

// Variant Calling and Consensus Generation
include { BWAMEM2_INDEX           } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM             } from '../../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_INDEX          } from '../../../modules/nf-core/samtools/index/main'
include { FREEBAYES               } from '../../../modules/local/freebayes/main'
include { PROCESS_GVCF            } from '../../../modules/local/process_gvcf/main'
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
    ch_reference           //
    ch_input_fastqs        // 
    ch_primer_bed          // value: path to reference file or null
    

    main:
    ch_versions = Channel.empty()

    //
    FASTP(
        ch_input_fastqs,
        [],
        '',
        '',
        ''
    )

    //
    BWAMEM2_INDEX(
        ch_reference
    )

    BWAMEM2_MEM(
        ch_input_fastqs,
        BWAMEM2_INDEX.out.index,
        ch_reference,
        'sort'
    )

    SAMTOOLS_INDEX(
        BWAMEM2_MEM.out.bam
    )
    ch_bam_bai = BWAMEM2_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0])

    if( ch_primer_bed ) {
        // IVAR Trim
    }

    FREEBAYES(
        ch_bam_bai,
        ch_reference
    )

    PROCESS_GVCF(
        FREEBAYES.out.gvcf,
        ch_reference
    )


    // For some reason it did not want to let me join the [] alone?
    //  So we've done it this way and that works    
    PROCESS_GVCF.out.ambiguous_vcf
        .combine(ch_reference.map{it -> it[1]})
        .map{ it -> [it[0], it[1], it[2], it[3], []]}
        .set{ ch_ambiguous_vcf_restructured }


    BCFTOOLS_CONSENSUS_AMBIGUOUS(
        ch_ambiguous_vcf_restructured
    )

    BCFTOOLS_CONSENSUS_FINAL(
        PROCESS_GVCF.out.consensus_vcf
            .join(BCFTOOLS_CONSENSUS_AMBIGUOUS.out.fasta, by: [0])
            .join(PROCESS_GVCF.out.mask, by: [0])
    )

    ADJUST_FASTA_HEADER(
        BCFTOOLS_CONSENSUS_FINAL.out.fasta,
        ch_reference,
        '.consensus',
        ''
    )

    emit:
    bam_bai   = ch_bam_bai
    consensus = ADJUST_FASTA_HEADER.out.consensus
    vcf       = PROCESS_GVCF.out.consensus_vcf
    versions  = ch_versions
}
