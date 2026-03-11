#!/usr/bin/env python3
'''Create Sample QC CSV file based on all pipeline outputs'''

import argparse
import csv
import json
import re
import statistics
import subprocess
import vcf

import pandas as pd
from Bio import SeqIO, SeqRecord
from pathlib import Path
from typing import Tuple


def init_parser() -> argparse.ArgumentParser:
    """
    Specify command line arguments
    Returns command line parser with inputs
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-b',
        '--bam',
        required=True,
        type=Path,
        help='Input sample bam file'
    )
    parser.add_argument(
        '-c1',
        '--consensus',
        required=True,
        type=Path,
        help='Input sample consensus sequence fasta file'
    )
    parser.add_argument(
        '-c2',
        '--consensus_n450',
        required=True,
        type=Path,
        help='Input sample N450 consensus sequence fasta file'
    )
    parser.add_argument(
        '-d',
        '--depth',
        required=True,
        type=Path,
        help='Input sample depth bed file'
    )
    parser.add_argument(
        '-n1',
        '--nextclade_n450',
        required=True,
        type=Path,
        help='Nextclade CSV file from Measles N450 dataset'
    )
    parser.add_argument(
        '-n2',
        '--nextclade_custom',
        required=True,
        type=Path,
        help='Nextclade CSV file from Measles Custom dataset'
    )
    parser.add_argument(
        '-s',
        '--sample',
        required=True,
        type=str,
        help='Sample name'
    )
    parser.add_argument(
        '-i',
        '--irida_id',
        required=True,
        type=str,
        help='IRIDA sample name to add to output'
    )
    parser.add_argument(
        '-g',
        '--genotype',
        required=True,
        type=str,
        help='Refernce genotype'
    )
    parser.add_argument(
        '-v',
        '--vcf',
        required=True,
        type=Path,
        help='Input sample passing vcf file'
    )
    parser.add_argument(
        '-m',
        '--matched_dsid',
        required=True,
        type=Path,
        help='TSV containing compared DSID calls'
    )
    parser.add_argument(
        '--nanoq_json',
        required=False,
        type=Path,
        help='Nanoq stats json file'
    )
    parser.add_argument(
        '--fastp_json',
        required=False,
        type=Path,
        help='Fastp stats json file'
    )
    parser.add_argument(
        '--seq_bed',
        required=False,
        type=Path,
        help='Input Sequencing primer bed file to test variants against'
    )
    parser.add_argument(
        '-r',
        '--ref_id',
        required=True,
        type=str,
        help='Name of Reference used'
    )
    return parser


def get_read_count(bam: Path) -> int:
    '''
    Purpose:
    --------
    Get the number of aligned reads from given bamfile using samtools and subprocess

    Parameters:
    -----------
    bam - Path
        Path to the primertrimmed bam file

    Returns:
    --------
    Integer number of aligned reads
    '''
    cmd = ['samtools', 'view', '-c', '-F0x900', bam]
    read_count = subprocess.run(cmd, capture_output=True, check=True, text=True).stdout.strip('\n')
    return int(read_count)


def parse_nanoq_json(injson: Path) -> float:
    '''
    Purpose:
    --------
    Parse nanoq json file to get wanted stats. Currently just the number of reads kept

    Parameters:
    -----------
    injson - Path
        Path to input nanoq json file

    Returns:
    --------
    Integer number of input reads
    '''
    with open(injson, 'r') as handle:
        json_obj = json.load(handle)
        return json_obj['reads']


def parse_fastp_json(injson: Path) -> float:
    '''
    Purpose:
    --------
    Parse fastp json file to get wanted stats. Currently just the number of reads kept

    Parameters:
    -----------
    injson - Path
        Path to input fastp json file

    Returns:
    --------
    Integer number of input reads
    '''
    with open(injson, 'r') as handle:
        json_obj = json.load(handle)
        return json_obj['filtering_result']['passed_filter_reads']


def parse_depth_bed(bed: Path, segment: range) -> Tuple[float, float, float]:
    '''
    Purpose:
    --------
    Calculate the mean and median sequencing depth from samtools depth bed file

    Parameters:
    -----------
    bed - Path
        Path to the input depth bed file from samtools depth

    Returns:
    --------
    Float mean sequencing depth
    Float median sequencing depth
    '''
    depth = []
    with open(bed, 'r') as handle:
        reader = csv.reader(handle, delimiter='\t')
        # Format is: Chrom, Pos, Depth
        above_20_count = 0
        for row in reader:
            if int(row[1]) in segment:
                depth.append(int(row[2]))

            if int(row[2]) >= 20 and int(row[1]) in segment:
                above_20_count += 1

    # Empty file return nothing
    if depth == []:
        return 0, 0, 0
    # Otherwise calc
    mean_dep = round(statistics.mean(depth), 2)
    median_dep = round(statistics.median(depth), 1)
    above_20_complete = round((above_20_count / len(depth) * 100), 2)
    return mean_dep, median_dep, above_20_complete


def parse_consensus(fasta: SeqRecord) -> Tuple[int, float, int, str]:
    '''
    Purpose:
    --------
    Parse consensus file to get the metrics based on it

    Parameters:
    -----------
    fasta - SeqRecord
        Bio SeqRecord of input consensus sequence

    Returns:
    --------
    List of 1-index N positions
    Integer number of Ns
    Float genome completeness calculated from the number of Ns
    Integer genome length
    String PASS if divisible, FAIL if not
    '''
    n_pos =  [i for i, base in enumerate(fasta.seq.lower(), start=1) if base == 'n']
    count_n = len(n_pos)
    seq_len = len(fasta.seq)
    completeness = (1 - (count_n / seq_len)) * 100
    completeness = round(completeness, 2)
    divisible = "PASS"
    if seq_len % 6 != 0:
        divisible = "FAIL"
    return n_pos, count_n, completeness, seq_len, divisible


def _create_variantpos_dict(var: Path, var_range: range) -> dict:
    '''Create dict with keys "variant" and "range"'''
    return {
        'variant': var,
        'range': var_range
    }


def parse_vcf(vcf_file: Path, n_pos: list) -> Tuple[str, list, str, dict]:
    '''
    Purpose:
    --------
    Parse input VCF file to find variants and their locations

    Parameters:
    -----------
    vcf_file - Path
        Path to input gzipped vcf file from args
    n_pos - List
        List of N positions to remove any low-depth filtered variants from reporting

    Returns:
    --------
    String of parsed variants to report
    List of variant-range dicts
    Dict of different variant counts
        {
            'total_variants': int,
            'num_snps': int,
            'num_iupac': int,
            'num_deletions': int,
            'num_deletion_sites': int,
            'num_insertions': int,
            'num_insertion_sites': int
        }
    '''
    # Base outputs
    variants = []
    variant_positions = []
    var_count_dict = {
        'total_variants': 0,
        'num_snps': 0,
        'num_iupac': 0,
        'num_deletions': 0,
        'num_deletion_sites': 0,
        'num_insertions': 0,
        'num_insertion_sites': 0
    }

    # Open and handle file
    with open(vcf_file, 'rb') as handle:
        reader = vcf.Reader(handle)
        for record in reader:
            # Odd issue previously - skip over Ns in vcf
            if str(record.ALT[0]).upper() == 'N':
                continue
            # Other odd issue, skip positions where alt is None
            if record.ALT[0] is None:
                continue
            # Make sure the variant wasn't masked by depth mask
            if int(record.POS) in n_pos:
                continue
            # Multiple alleles not supported warning
            #  These should have been corrected in previous processing steps as well to only have 1
            if len(record.ALT) > 1:
                print(f'WARNING: Multiple alleles not supported currently. Using only the first one for position: {record.POS}')

            # Create string of variant and add to list of variants along with getting the lengths of the ref and alt
            ref_len = len(record.REF)
            alt_len = len(record.ALT[0])
            alt_str = str(record.ALT[0])
            # For IUPACs make sure that the info column exists and is ambiguous
            #  This is just for Illumina, no Nanopore ambiguous at the moment
            iupac = False
            if ("ConsensusTag" in record.INFO) and (record.INFO["ConsensusTag"] == 'ambiguous'):
                iupac = True
                alt_str = record.INFO["ConsensusBase"]
            variant = f'{record.REF}{record.POS}{alt_str}'

            # Type of mutation leads to different spots affected and tracking
            if variant not in variants:
                # Deletions
                if ref_len > alt_len:
                    # For dels, POS is the kept genomic position so +1 to start to get the deleted positions
                    # The range works as a 9bp deletion would be length of 10 in vcf record
                    del_range = range(record.POS+1,record.POS+ref_len)
                    variants.append(variant)
                    variant_positions.append(_create_variantpos_dict(variant, del_range))
                    var_count_dict['num_deletions'] += len(del_range)
                    var_count_dict['num_deletion_sites'] += 1

                # Insertions
                elif ref_len < alt_len:
                    variants.append(variant)
                    variant_positions.append(_create_variantpos_dict(variant, range(record.POS, record.POS+1)))
                    var_count_dict['num_insertions'] += (alt_len - ref_len) # As ref length can differ
                    var_count_dict['num_insertion_sites'] += 1

                # Multiple SNPs together
                elif (ref_len > 1) and (ref_len == alt_len):
                    mult_snp_range = range(record.POS, record.POS+len(record.REF))
                    for i, ref_base in enumerate(record.REF):
                        # Include only snps incase the variant is like ref=ATG alt=TTC where the T isn't a SNP
                        if ref_base == alt_str[i]:
                            continue
                        variant = f'{ref_base}{mult_snp_range[i]}{alt_str[i]}'
                        variants.append(variant)
                        variant_positions.append(_create_variantpos_dict(variant, range(mult_snp_range[i], mult_snp_range[i]+1)))
                        if iupac:
                            var_count_dict['num_iupac'] += 1
                        else:
                            var_count_dict['num_snps'] += 1
                else:
                    variants.append(variant)
                    variant_positions.append(_create_variantpos_dict(variant, range(record.POS, record.POS+1)))
                    if iupac:
                        var_count_dict['num_iupac'] += 1
                    else:
                        var_count_dict['num_snps'] += 1

    # Final Summary and Return
    if variants:
        var_count_dict['total_variants'] = len(variants)
        return ';'.join(variants), variant_positions, var_count_dict
    return 'none', variant_positions, var_count_dict


def range_overlap(r1: range, r2: range) -> bool:
    '''Return True if range2 overlaps range1'''
    x1, x2 = r1.start, r1.stop
    y1, y2 = r2.start, r2.stop
    return x1 <= y2 and y1 <= x2


def check_primers(bed: Path, variant_locations: list) -> str:
    '''
    Purpose:
    --------
    Parse input bed file for any variant regions that overlap

    Parameters:
    -----------
    bed - Path
        Path to bed file
    variant_locations - list[dict]
        List of variantpos dictionaries with variant and range keys

    Returns:
    --------
    String of primer mutations found or string "none" if there are none
    '''
    if not variant_locations:
        return 'none'

    primer_mutations = []
    with open(bed, 'r') as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            # Bed file needs at least 4 rows (chrom, start, stop, name)
            if len(row) < 4:
                continue
            # Set primer values, make sure start lower than stop for range
            start, stop, name = int(row[1]), int(row[2]), row[3]
            if start > stop:
                start, stop = stop, start
            location = range(start, stop + 1) # Plus one to make sure that we get mutations in the final location of the range

            # Check if the location range overlaps with any variant ranges
            for var_dict in variant_locations:
                if range_overlap(location, var_dict['range']):
                    primer_mutations.append(f'{var_dict["variant"]}-{name}')

    if not primer_mutations:
        return 'none'
    return ';'.join(primer_mutations)


def _get_nextclade_row_dict(nextclade_csv: Path, sample: str) -> dict:
    '''
    Purpose:
    --------
    Parse nextclade CSV file to return row-dict

    Parameters:
    -----------
    nextclade_csv - Path
        Path to nextclade CSV file. ';' delimited
    sample - str
        String sample name to match to nextclade CSV

    Returns:
    --------
    Dict values of wanted nextclade row
    '''
    with open(nextclade_csv, 'r') as handle:
        reader = csv.DictReader(handle, delimiter=';')
        # There should only be 1 sample / 1 line but just in case
        for d in reader:
            seqname = d['seqName'].split(' ')[0]
            if seqname == sample:
                return d
    return {}


def parse_n450_nextclade(nextclade_csv: Path, sample: str) -> Tuple[str, range]:
    '''
    Purpose:
    --------
    Parse nextclade N450 CSV file for genotype and other N450 metrics

    Parameters:
    -----------
    nextclade_csv - Path
        Path to nextclade CSV file. ';' delimited
    sample - str
        String sample name to match to nextclade CSV

    Returns:
    --------
    String nextclade clade
    Range where N450 is located
    '''
    d = _get_nextclade_row_dict(nextclade_csv, sample)
    if d:
        genotype = d['clade']
        # insertions structured as 0:BASES,450:BASES as all alignments will be 450 bp with the dataset used
        to_n450 = 0
        if ',' in d['insertions']:
            to_n450 = len(d['insertions'].split(',')[0].split(':')[1])
        return genotype, range(to_n450, to_n450+451)
    else:
        return '', range(0,452)


def get_custom_nextclade_vals(nextclade_csv: Path, sample: str) -> Tuple[str, str, str]:
    '''
    Purpose:
    --------
    Parse custom nextclade CSV file to information on potential issue sites

    Parameters:
    -----------
    nextclade_csv - Path
        Path to nextclade CSV file. ';' delimited
    sample - str
        String sample name to match to nextclade CSV

    Returns:
    --------
    Tuple of 3 strings identifying frameshifts, stop codons, and mutated stop codons
    '''
    # Nextclade CDS Checks
    d = _get_nextclade_row_dict(nextclade_csv, sample)
    if d:
        aa_mutations = d['aaSubstitutions']
        frameshifts = d['qc.frameShifts.frameShifts']
        stop_codons = d['qc.stopCodons.stopCodons']

        stop_codon_pattern = '|'.join(["N:\\*", "C:\\*", "P:\\*", "V:\\*", "M:\\*", "F:\\*", "H:\\*", "L:\\*"])
        mutated_stop_codons_match = re.findall(stop_codon_pattern, aa_mutations)
        mutated_stop_codons = '|'.join(mutated_stop_codons_match)

        return frameshifts, stop_codons, mutated_stop_codons
    else:
        return '', '', ''


def parse_dsid_tsv(matched_dsid: Path, sample: str) -> str:
    '''
    Purpose:
    --------
    Parse matched dsid tsv file to determine if the sample matched a known dsid and get N450 completeness

    Parameters:
    -----------
    matched_dsid - Path
        Path to matched dsid TSV
    sample - str
        String sample name to match to dsid file

    Returns:
    --------
    String of DSID or its value in input table or No Data
    Float N450 completeness
    '''
    # Check for a match
    with open(matched_dsid, 'r') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for d in reader:
            if d['sample'] == sample:
                match = str(d['matched_dsid'])
                n450_completeness = float(d['completeness'])
                return match, n450_completeness
    return 'No Data', 0.00


def grade_n450(dsid: str, n450_mean_depth: float, n450_completeness: float) -> str:
    '''
    Purpose:
    --------
    Determine if the sample N450 passes internal QC metrics

    Parameters:
    -----------
    dsid - str
        DSId value determined
    n450_mean_depth - float
        Mean n450 sequencing depth
    n450_completeness - float
        How complete the N450 region is

    Returns:
    --------
    String N450 QC status of PASS / FAIL
    '''
    if (dsid in ['No Data', 'Incomplete', 'Ambiguous Base']) or (n450_mean_depth < 30) or (n450_completeness < 100):
        return 'FAIL'
    return 'PASS'


def grade_qc(completeness: float, mean_dep: float, median_dep: float, divisible: str,
             frameshift: bool, nonsense_mutation: bool, stop_mutation: bool, genotype_match: bool) -> str:
    '''
    Purpose:
    --------
    Determine if the sample passes internal QC metrics

    Parameters:
    -----------
    completeness - float
        Genome completeness
    mean_dep - float
        Mean sequencing depth
    median_dep - float
        Median sequencing depth
    divisible - str
        If sample is divisible by 6 = PASS
    frameshift - bool
        If there was a frameshift mutation identified
    nonsense_mutation - bool
        If there was a nonsense mutation identified
    stop_mutation - bool
        If there was a stop codon mutated
    genotype_match - bool
        If the reference genotype matches the samples
    nextclade_csv - Path
        Path to custom WG nextclade CSV to determine metrics from

    Returns:
    --------
    String QC status
    '''
    qc_status = []
    # Completeness
    if completeness < 90:
        if completeness < 50:
            return 'INCOMPLETE_GENOME'
        else:
            qc_status.append('PARTIAL_GENOME')
    # Overall Coverage Depth
    if (mean_dep < 20) or (median_dep < 20):
        qc_status.append('LOW_SEQ_DEPTH')
    # Divisible by 6
    if divisible == "FAIL":
        qc_status.append('NOT_DIVISIBLE')
    # Frameshift
    if frameshift:
        qc_status.append('FRAMESHIFT_MUTATION')
    # Nonsense
    if nonsense_mutation:
        qc_status.append('NONSENSE_MUTATION')
    # Stop
    if stop_mutation:
        qc_status.append('STOP_CODON_MUTATION')
    # Reference genotype match
    if not genotype_match:
        qc_status.append('UNMATCHING_GENOTYPE')

    if qc_status:
        return ';'.join(qc_status)
    return 'PASS'


def main() -> None:
    '''Run the program'''
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Final read count information based on pipeline filtering tool
    if args.nanoq_json:
        num_input_reads = parse_nanoq_json(args.nanoq_json)
    elif args.fastp_json:
        num_input_reads = parse_fastp_json(args.fastp_json)
    else:
        num_input_reads = 0

    num_aligned_reads = get_read_count(args.bam)
    consensus = SeqIO.read(args.consensus, "fasta")
    n_pos, count_n, completeness, seq_len, divisible = parse_consensus(consensus)
    mean_dep, median_dep, above_20_completeness = parse_depth_bed(args.depth, range(0,seq_len+1)) # Use the whole genome for calc
    variants, variant_positions, var_count_dict = parse_vcf(args.vcf, n_pos)
    frameshift, nonsense, stop_mutation = get_custom_nextclade_vals(args.nextclade_custom, args.sample)

    # Optional inputs

    seq_primer_overlap = 'NA'
    if args.seq_bed:
        seq_primer_overlap = check_primers(args.seq_bed, variant_positions)

    # N450 dataset and DSId checks
    matched_dsid, n450_completeness = parse_dsid_tsv(args.matched_dsid, args.sample)

    genotype, n450_range = parse_n450_nextclade(args.nextclade_n450, args.sample)
    n450_mean_depth = 0
    n450_median_depth = 0

    if n450_completeness > 0:
        n450_mean_depth, n450_median_depth, _ = parse_depth_bed(args.depth, n450_range)
    n450_status = grade_n450(matched_dsid, n450_mean_depth, n450_completeness)

    # MF-NCR
    mf_mean, mf_med, mf_20x = parse_depth_bed(args.depth, range(4338,5350))

    # Grade qc
    frameshift_status = (frameshift != '')
    nonsense_status = (nonsense != '')
    stop_mutation_status = (stop_mutation != '')
    genotype_match = (args.genotype == genotype)
    qc_status = grade_qc(completeness, mean_dep, median_dep, divisible,
                         frameshift_status, nonsense_status, stop_mutation_status,
                         genotype_match)

    # N450 seq for inclusion in final excel output ONLY
    #  If empty file, get value error and then put in 450 Ns
    try:
        n450_seq = SeqIO.read(args.consensus_n450, "fasta").seq.upper()
    except ValueError:
        n450_seq = 'N'*450

    # Output
    final = {
        'sample': [args.sample],
        'genotype': [genotype],
        'reference': [args.ref_id],
        'matched_dsid': [matched_dsid],
        'num_input_reads': [num_input_reads],
        'num_aligned_reads': [num_aligned_reads],
        'num_consensus_n': [count_n],
        'genome_completeness_percent': [completeness],
        'Above 20x Coverage': [above_20_completeness],
        'mean_sequencing_depth': [mean_dep],
        'median_sequencing_depth': [median_dep],
        'total_variants': [var_count_dict['total_variants']],
        'num_snps': [var_count_dict['num_snps']],
        'num_iupac': [var_count_dict['num_iupac']],
        'num_deletions': [var_count_dict['num_deletions']],
        'num_deletion_sites': [var_count_dict['num_deletion_sites']],
        'num_insertions': [var_count_dict['num_insertions']],
        'num_insertion_sites': [var_count_dict['num_insertion_sites']],
        'genome_length': [seq_len],
        'divisible_by_6': [divisible],
        'frameshifts': [frameshift],
        'nonsense_mutation': [nonsense],
        'mutated_stop_codon': [stop_mutation],
        'variants': [variants],
        'sequencing_primer_variants': [seq_primer_overlap],
        'N450_completeness': [n450_completeness],
        'N450_mean_depth': [n450_mean_depth],
        'N450_status': [n450_status],
        'MF-NCR Mean': [mf_mean],
        'MF-NCR Median': [mf_med],
        'qc_status': [qc_status],
        'N450_fasta': [f">{args.sample}-N450\n{n450_seq}"],
        'genome_fasta': [f">{args.sample}\n{consensus.seq.upper()}"],
        'irida_id': [args.irida_id]
    }
    df = pd.DataFrame.from_dict(final)
    df.to_csv(f'{args.sample}.qc.csv', index=False)

if __name__ == '__main__':
    main()
