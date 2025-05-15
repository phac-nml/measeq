#!/usr/bin/env python3
'''Create Sample QC CSV file based on pipeline outputs'''

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
        '-c',
        '--consensus',
        required=True,
        type=Path,
        help='Input sample consensus sequence file'
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
        '-r',
        '--read_json',
        required=False,
        type=Path,
        help='Nanoq stats json file'
    )
    parser.add_argument(
        '-s',
        '--sample',
        required=True,
        type=str,
        help='Sample name'
    )
    parser.add_argument(
        '-s2',
        '--strain',
        required=True,
        type=str,
        help='Refernce Stain'
    )
    parser.add_argument(
        '-v',
        '--vcf',
        required=True,
        type=Path,
        help='Input sample passing vcf file'
    )
    parser.add_argument(
        '--seq_bed',
        required=False,
        type=Path,
        help='Input Sequencing primer bed file to test variants against'
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

def parse_depth_bed(bed: Path) -> Tuple[float, float]:
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
        for row in reader:
            depth.append(int(row[2]))

    # Empty file return nothing
    if depth == []:
        return 0, 0
    # Otherwise calc
    mean_dep = round(statistics.mean(depth), 2)
    median_dep = round(statistics.median(depth), 1)
    return mean_dep, median_dep

def parse_consensus(fasta: SeqRecord) -> Tuple[int, float, int, bool]:
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
    Integer number of Ns
    Float genome completeness calculated from the number of Ns
    Integer genome length
    Bool TRUE if divisible, FALSE if not
    '''
    n_pos =  [i for i, base in enumerate(fasta.seq.lower()) if base == 'n']
    count_n = len(n_pos)
    seq_len = len(fasta.seq)
    completeness = (1 - (count_n / seq_len)) * 100
    completeness = round(completeness, 2)
    divisible = True
    if seq_len % 6 != 0:
        divisible = False
    return count_n, completeness, seq_len, divisible

def _create_variantpos_dict(var: Path, var_range: range) -> dict:
    '''Create dict with keys "variant" and "range"'''
    return {
        'variant': var,
        'range': var_range
    }

def parse_vcf(vcf_file: Path) -> Tuple[str, list, str, dict]:
    '''
    Purpose:
    --------
    Parse input VCF file to find variants and their locations

    Parameters:
    -----------
    vcf_file - Path
        Path to input gzipped vcf file from args

    Returns:
    --------
    String of parsed variants to report
    List of variant-range dicts
    Dict of different variant counts
        {'total_variants': int, 'num_snps': int, 'num_deletions': int, 'num_deletion_sites': int, 'num_insertions': int, 'num_insertion_sites': int}
    '''
    # Base outputs
    variants = []
    variant_positions = []
    var_count_dict = {
        'total_variants': 0,
        'num_snps': 0,
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
            if record.ALT[0] == None:
                continue
            # Multiple alleles not supported warning
            if len(record.ALT) > 1:
                print(f'WARNING: Multiple alleles not supported currently. Using only the first one for position: {record.POS}')

            # Create string of variant and add to list of variants along with getting the lengths of the ref and alt
            variant = f'{record.REF}{record.POS}{record.ALT[0]}'
            ref_len = len(record.REF)
            alt_len = len(record.ALT[0])
            alt_str = str(record.ALT[0])

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
                    var_count_dict['num_insertions'] += (alt_len - 1) # -1 for the included ref base
                    var_count_dict['num_insertion_sites'] += 1
                # Multiple SNPs together
                elif (ref_len > 1) and (ref_len == alt_len):
                    mult_snp_range = range(record.POS, record.POS+len(record.REF))
                    for i, ref_base in enumerate(record.REF):
                        # Include only snps incase the variant is like ref=ATG alt=TTC where the T isn't a SNP
                        if ref_base == alt_str[i]:
                            pass
                        variant = f'{ref_base}{mult_snp_range[i]}{alt_str[i]}'
                        variants.append(variant)
                        variant_positions.append(_create_variantpos_dict(variant, range(mult_snp_range[i], mult_snp_range[i]+1)))
                        var_count_dict['num_snps'] += 1
                else:
                    variants.append(variant)
                    variant_positions.append(_create_variantpos_dict(variant, range(record.POS, record.POS+1)))
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

def get_strain(nextclade_csv: Path, sample: str) -> str:
    '''
    Purpose:
    --------
    Parse nextclade N450 CSV file to get the strain

    Parameters:
    -----------
    nextclade_csv - Path
        Path to nextclade CSV file. ';' delimited
    sample - str
        String sample name to match to nextclade CSV

    Returns:
    --------
    String nextclade clade
    '''
    d = _get_nextclade_row_dict(nextclade_csv, sample)
    if d:
        return d['clade']
    else:
        return ''

def get_custom_nextclade_vals(nextclade_csv: Path, sample: str) -> Tuple[str, str, str]:
    '''
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

def grade_qc(completeness: float, mean_dep: float, median_dep: float, divisible: bool, 
             frameshift: bool, nonsense_mutation: bool, stop_mutation: bool, strain_match: bool) -> str:
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
    divisible - bool
        If sample is divisible by 6
    frameshift - bool
        If there was a frameshift mutation identified
    nonsense_mutation - bool
        If there was a nonsense mutation identified
    stop_mutation - bool
        If there was a stop codon mutated
    strain_match - bool
        If the reference strain matches the samples
    nextclade_csv - Path
        Path to custom WG nextclade CSV to determine metrics from

    Returns:
    --------
    String QC status
    '''
    qc_status = []
    # Completeness
    if completeness < 0.9:
        if completeness < 0.5:
            return 'INCOMPLETE_GENOME'
        else:
            qc_status.append('PARTIAL_GENOME')
    # Coverage Depth
    if (mean_dep < 20) or (median_dep < 20):
        qc_status.append('LOW_SEQ_DEPTH')
    # Divisible by 6
    if not divisible:
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
    # Reference strain match
    if not strain_match:
        qc_status.append('UNMATCHING_STRAIN')

    if qc_status:
        return ';'.join(qc_status)
    return 'PASS'

def main() -> None:
    '''Run the program'''
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Do something with the data
    strain = get_strain(args.nextclade_n450, args.sample)
    #num_input_reads = parse_nanoq_json(args.read_json)
    num_input_reads = 1000
    num_aligned_reads = get_read_count(args.bam)
    consensus = SeqIO.read(args.consensus, "fasta")
    count_n, completeness, seq_len, divisible = parse_consensus(consensus)
    mean_dep, median_dep = parse_depth_bed(args.depth)
    variants, variant_positions, var_count_dict = parse_vcf(args.vcf)
    frameshift, nonsense, stop_mutation = get_custom_nextclade_vals(args.nextclade_custom, args.sample)

    # Optional inputs
    seq_primer_overlap = 'NA'
    if args.seq_bed:
        seq_primer_overlap = check_primers(args.seq_bed, variant_positions)

    # Grade qc
    frameshift_status = (frameshift != '')
    nonsense_status = (nonsense != '')
    stop_mutation_status = (stop_mutation != '')
    strain_match = (args.strain == strain)
    qc_status = grade_qc(completeness, mean_dep, median_dep, divisible,
                         frameshift_status, nonsense_status, stop_mutation_status,
                         strain_match)

    # Output
    final = {
        'sample': [args.sample],
        'strain': [strain],
        'num_input_reads': [num_input_reads],
        'num_aligned_reads': [num_aligned_reads],
        'num_consensus_n': [count_n],
        'genome_completeness_percent': [completeness],
        'mean_sequencing_depth': [mean_dep],
        'median_sequencing_depth': [median_dep],
        'total_variants': [var_count_dict['total_variants']],
        'num_snps': [var_count_dict['num_snps']],
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
        'qc_status': [qc_status]
    }
    df = pd.DataFrame.from_dict(final)
    df.to_csv(f'{args.sample}.qc.csv', index=False)

if __name__ == '__main__':
    main()
