#!/usr/bin/env python
'''Script to create amplicon bed file from primer bed file with minimal dependencies'''

import argparse
import csv
import re
from collections import defaultdict
from typing import Generator

def init_parser() -> argparse.ArgumentParser:
    """
    Specify command line arguments
    Returns command line parser with inputs
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-b',
        '--bed',
        required=True,
        type=str,
        help='Path to input primer bed file'
    )
    parser.add_argument(
        '-o',
        '--outfile',
        required=False,
        default='amplicon.bed',
        type=str,
        help='Name for output amplicon bed file. Default: "amplicon.bed"'
    )
    return parser

def primer_pair_generator(bed: str, primer_rgx=r'^(.*)_(LEFT|RIGHT).*$') -> Generator[str, list, list]:
    '''
    Generate and yield primer pairs based on the primer name
        Primers are paired through the primer_rgx which requires the names to contain:
            _LEFT or _RIGHT
    '''
    pair_dict = defaultdict(lambda: [None, None])
    with open(bed, 'r') as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            # Bed file needs at least 5 rows (chrom, start, stop, name, pool)
            if len(row) < 5:
                continue
            name = str(row[3])

            # Pair based on if LEFT or RIGHT match and yield upon a pair
            primer_match = re.search(primer_rgx, name)
            if primer_match:
                if primer_match.group(2) == 'LEFT':
                    primer_prefix = primer_match.group(1)
                    if primer_prefix not in pair_dict:
                        pair_dict[primer_prefix][0] = row
                    # In case of alt primers, don't want to yield yet
                    elif pair_dict[primer_prefix][0] != None:
                        print(f'Skipping {name} - already has a primer in dict')
                        continue
                    else:
                        yield primer_prefix, pair_dict[primer_prefix][0], row
                        del pair_dict[primer_prefix]
                elif primer_match.group(2) == 'RIGHT':
                    primer_prefix = primer_match.group(1)
                    if primer_prefix not in pair_dict:
                        pair_dict[primer_prefix][1] = row
                    # In case of alt primers, don't want to yield yet
                    elif pair_dict[primer_prefix][1] != None:
                        print(f'Skipping {name} - already has a primer in dict')
                        continue
                    else:
                        yield primer_prefix, pair_dict[primer_prefix][0], row
                        del pair_dict[primer_prefix]


def main() -> None:
    '''Run the program'''
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()
    # full_genome_region = defaultdict(lambda: defaultdict(str))
    # Use variables to find the starting and end position of the scheme
    scheme_start = 1000000000
    scheme_end = 0

    # Write amplicon bed file
    with open(args.outfile, 'w') as f:
        for ppair in primer_pair_generator(args.bed):
            # ppair = [ primer_name, [primer-left] , [primer-right] ]
            # Check starts and stops to make sure that start < stop and if not switches them
            start_l, stop_l, start_r, stop_r = int(ppair[1][1]), int(ppair[1][2]), int(ppair[2][1]), int(ppair[2][2])
            if start_l > stop_l:
                start_l, stop_l = stop_l, start_l
            if start_r > stop_r:
                start_r, stop_r = stop_r, start_r

            # Track overall genome start and end
            if scheme_start > stop_l:
                scheme_start = stop_l
            if scheme_end < start_r:
                scheme_end = start_r

            # Write output
            chrom = ppair[1][0]
            pname = ppair[0]
            pool = ppair[1][4]
            f.write(f'{chrom}\t{stop_l}\t{start_r}\t{pname}\t{pool}\t+\n')

    # Write the region tiled to file
    with open(f'tiling_region.bed', 'w') as f:
        f.write(f'{chrom}\t{scheme_start}\t{scheme_end}\n')

if __name__ == '__main__':
    main()
