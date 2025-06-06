#!/usr/bin/env python3

import argparse
import pysam
import pandas as pd
import re
import sys
from pathlib import Path

from collections import Counter

def init_parser() -> argparse.ArgumentParser:
    '''
    Purpose
    -------
    Parse CL inputs to be used in script
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-b',
        '--bam',
        required=True,
        help='Path to input BAM file. Make sure that the BAI file is also on the same path'
    )
    parser.add_argument(
        '-r',
        '--reference',
        required=True,
        help='Input reference fasta file'
    )
    parser.add_argument(
        '--sample',
        required=False,
        help='Input sample name to be used in output file creation'
    )
    parser.add_argument(
        '--min_base_qual',
        required=False,
        type=int,
        default=10,
        help='Minimum base quality to use read for position'
    )
    parser.add_argument(
        '--min_map_qual',
        required=False,
        type=int,
        default=30,
        help='Minimum mapping quality to use read for position'
    )
    parser.add_argument(
        '--min_read_count',
        required=False,
        type=int,
        default=10,
        help='Minimum number of reads to report position'
    )
    parser.add_argument(
        '--min_report_percent',
        required=False,
        type=int,
        default=15,
        help='Minimum percentage of non-reference bases to report out'
    )
    parser.add_argument(
        '--variant_call_threshold',
        required=False,
        type=int,
        default=75,
        help='Minimum percentage for a variant to called'
    )

    return parser

def create_ref_dict(ref_fasta: str) -> dict:
    """
    Purpose
    -------
    Process reference fasta file with pysam to create a dictionary of indexs with their bases
        It is 0-based

    Parameters
    ----------
    ref_fasta : str
        Path to input reference fasta file

    Returns
    -------
    dict{ Position: Base }
    """
    ref_dict = {}
    ref_fasta = pysam.FastaFile(ref_fasta)
    for idx, base in enumerate(ref_fasta.fetch(ref_fasta.references[0])):
        # Allows lowercase and uppercase refs, fixes small bug
        base = base.upper()
        ref_dict[idx] = base

    return ref_dict

def ensure_all_bases_dict_from_counter(counter: Counter) -> dict:
    """Write a summary"""
    d = {base: counter.get(base, 0) for base in ['A', 'T', 'G', 'C', '*', 'ins']}
    d['del'] = d.pop('*')
    return d

def calculate_non_ref_percent(pileup_dict: dict, ref_base: str) -> float:
    """
    Purpose
    -------
    Calculate the non-reference base percentage for the given position

    Parameters
    ----------
    pileup_dict: dict
        Contains read counts for the current position
    ref_base: str
        Reference base

    Returns
    -------
    Non-reference percentage as float
    """
    total = pileup_dict["total_reads"]
    non_ref_percent = 0

    if ref_base != "N":
        ref_base_count = pileup_dict[ref_base]
        if total == 0:
            non_ref_percent = 0
        elif ref_base_count == 0:
            non_ref_percent = 100
        else:
            non_ref_percent = round((100 - ((ref_base_count / total) * 100)), 2)

    return non_ref_percent

def calc_base_percent(pileup_dict: dict, base: str) -> float:
    """
    Purpose
    -------
    Calculate the base percentage for the position

    Parameters
    ----------
    pileup_dict: dict
        Contains read counts for the current position
    base: str
        Dict key of the base to use for percent calculation

    Returns
    -------
    Given base percentage as float
    """
    return round(((pileup_dict[base] / pileup_dict["total_reads"]) * 100), 2)

def determine_type(pileup_dict: dict, most_common_base: str, variant_call_threshold: int) -> str:
    """
    Purpose
    -------
    Determine if the location is a mutation, mixed site, or reference base

    Parameters
    ----------
    ref_dict: dict
        Contains: {Position: Base} for the reference genome
    most_common_base: str
        The most common base found at pileup location
    variant_call_threshold: int
        Threshold at which a variant is called. Default: 70

    Returns
    -------
    String location status
    """
    if most_common_base not in ['A', 'T', 'C', 'G', 'del']:
        return 'Insertion'

    if (100 - variant_call_threshold) >= pileup_dict['percentage_nonref']:
        return 'Ref'
    elif calc_base_percent(pileup_dict, most_common_base) >= variant_call_threshold:
        if most_common_base != 'del':
            return f'{most_common_base} SNP'
        else:
            return 'Deletion'
    else:
        return 'Mixed'

def parse_variation_from_bam(bamfile: str, ref_dict: dict, base_q: int,
                             map_q: int, min_read_count: int,
                             min_report_percent: int, variant_call_threshold: int
                            ) -> list:
    """
    Purpose
    -------
    Parse BAM file using pysam mpileup to find positions where non-ref base count is above
        the min_report_percent percent and min_read_count

    Parameters
    ----------
    bamfile : str
        Path to input BAM file
    ref_dict: dict
        Contains: {Position: Base} for the reference genome
    base_q: int
        Minimum base quality to use read. Default: 10
    map_q: int
        Minimum mapping quality to use read. Default: 30
    min_read_count: int
        Minimum reads required to report position. Default: 20
    min_report_percent: int
        Minimum non-reference base percentage to report out. Default: 15
    variant_call_threshold: int
        Threshold at which a variant is called. Default: 70

    Returns
    -------
    list of dictionaries containing positional variation information if variation is above min_report_percent
    """
    bamfile = pysam.AlignmentFile(bamfile, "rb" )
    variation_info = []

    # Based on pysam docs for iterating over each base and Piranha's variation analysis
    for pileupcolumn in bamfile.pileup(bamfile.references[0], min_base_quality=base_q, min_mapping_quality=map_q):
        # Set vars to make access more concise and easier
        num_reads = pileupcolumn.get_num_aligned()
        zero_idx_pos = pileupcolumn.reference_pos
        ref_base = ref_dict[zero_idx_pos]

        # Exit as early as possible if it is NOT a variant that is > min_report_percent
        #  Also stops 0 read locations from making it in
        if num_reads < min_read_count:
            continue

        # For each read in the specific position, count total base found
        data_list = pileupcolumn.get_query_sequences(add_indels=True)
        case_insensitive_pos_counter = Counter(map(str.upper, data_list))

        # Positions before deletions have "-" in them with the deletion length along with the next positions
        #  containing *s for the deletions. We are going to summarize deletions basd on the *s so need
        #  to skip the spots before while tracking the base in the current position
        #  Ex. ['a-1n', 'a-1n', 'A-1N', 'A-1N']
        if any([x for x in data_list if '-' in x]):
            case_insensitive_pos_counter = Counter({k: c for k, c in case_insensitive_pos_counter.items() if '-' not in k})
            base_list = [x[0].upper() for x in data_list if '-' in x]
            case_insensitive_pos_counter = case_insensitive_pos_counter + Counter(base_list)

        # Can skip spots where we have the reference base as the most common and it isn't above the minimum report percent
        most_common_base, read_count = case_insensitive_pos_counter.most_common(1)[0]
        if (most_common_base == ref_base) and (read_count/num_reads * 100 > 100 - min_report_percent):
            continue

        # Analyze insertions differently
        # Insertions are hard as they can be quite variable and we want a summary for them
        #  Ex. ['A+6AAAAAG', 'A+6AAAAAG', 'a+5aaagg', 'a+6aaaaag']
        if (any([x for x in data_list if '+' in x])):
            ins_list = [x for x in data_list if '+' in x]
            case_insensitive_pos_counter['ins'] = len(ins_list)

            # Do something with it to summarize and adjust this later
            ins_pattern = re.compile(r"\d+")
            ins_len_counts = Counter((f'{ins_pattern.findall(x)[0]}bp' for x in ins_list))
            ins_summary_str = ", ".join(f"{length}: {count}" for length, count in ins_len_counts.most_common())

        # Add missing bases to the counter, change * to del, and make dict
        pileup_dict = ensure_all_bases_dict_from_counter(case_insensitive_pos_counter)

        # Add values to dictionary
        #  chrom is second as multiqc uses first column for output and having it be the same leads to issues
        #  if we are using multiqc that is
        pileup_dict["position"] = zero_idx_pos + 1
        pileup_dict["chrom"] = pileupcolumn.reference_name
        pileup_dict["total_reads"] = num_reads
        pileup_dict["percentage_nonref"] = calculate_non_ref_percent(pileup_dict, ref_base)
        pileup_dict["ref_base"] = ref_base

        # Type of mutation
        if most_common_base == '*':
            most_common_base = 'del'
        pileup_dict["variant_type"] = determine_type(pileup_dict, most_common_base, variant_call_threshold)

        # If INS < variant_call_threshold we don't want to add insertion summary data
        if calc_base_percent(pileup_dict, 'ins') < (100 - variant_call_threshold):
            ins_summary_str = 'NA'
        pileup_dict["distribution"] = ins_summary_str

        # Add to report
        variation_info.append(pileup_dict)

    return variation_info

def main() -> None:
    '''
    Main entry point into script
    '''
    ## Init Parser and Set Arguments ##
    parser = init_parser()
    args = parser.parse_args()

    # Current test analysis
    path_bam = Path(args.bam)
    if args.sample:
        sample_name = args.sample
    else:
        sample_name = path_bam.name.split('.')[0]
    if path_bam.stat().st_size < 200:
        print(f"{args.bam} is too small a file. Skipping")
        sys.exit(0)

    ref_dict = create_ref_dict(args.reference)
    var_info = parse_variation_from_bam(path_bam, ref_dict, args.min_base_qual,
                                        args.min_map_qual, args.min_read_count,
                                        args.min_report_percent, args.variant_call_threshold)
    if var_info == []:
        sys.exit(0)
    df = pd.DataFrame.from_dict(var_info)
    # Reorder cols the lazy way
    df = df[['position', 'chrom', 'A', 'C', 'T', 'G', 'del', 'ins',
             'total_reads', 'percentage_nonref', 'ref_base',
             'variant_type', 'distribution'
            ]]
    df.to_csv(f'{sample_name}_variation.csv', index=False)

if __name__ == "__main__":
    main()
