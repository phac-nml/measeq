#!/usr/bin/env python3
'''
Check how complete an amplicon is based on bed file and consensus
'''
import argparse
import csv
from Bio import SeqIO, SeqRecord

def init_parser() -> argparse.ArgumentParser:
    """
    Specify command line arguments
    Returns command line parser with inputs
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-s',
        '--sample',
        required=True,
        type=str,
        help='Sample name'
    )
    parser.add_argument(
        '-a',
        '--amplicon_bed',
        required=True,
        type=str,
        help='Path to amplicon bed file'
    )
    parser.add_argument(
        '-c',
        '--consensus',
        required=True,
        type=str,
        help='Input sample consensus sequence file'
    )
    return parser

def parse_consensus(fasta: SeqRecord) -> str:
    '''
    Purpose:
    --------
    Parse consensus file to get the genome completeness and N count

    Parameters:
    -----------
    fasta - SeqRecord
        Bio SeqRecord of input consensus sequence

    Returns:
    --------
    Integer number of Ns
    Float genome completeness calculated from the number of Ns
    '''
    n_pos =  [i for i, base in enumerate(fasta.seq.lower()) if base == 'n']
    count_n = len(n_pos)
    completeness = 1 - (count_n / len(fasta.seq))
    completeness = round(completeness, 4)
    return count_n, completeness

def main() -> None:
    '''Run the program'''
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Read consensus file
    consensus = SeqIO.read(args.consensus, "fasta")
    consensus.seq = consensus.seq.upper()

    # Read in the amplicon bed file
    #  It shouldn't change format as it comes from the primers_to_amplicons script
    with open(args.amplicon_bed, 'r') as handle:
        reader = csv.reader(handle, delimiter='\t')
        out = {}
        for row in reader:
            # Remember that by primer conventions start is 0-based and stop is 1-based
            start, stop, name = int(row[1]), int(row[2]), str(row[3])
            amp_length = len(range(start, stop))
            n_count = consensus.seq[start:stop-1].count('N') # -1 to stop to get to 0-based for correct position

            # Calc
            if amp_length == 0:
                continue
            else:
                completeness = round(1 - (n_count/amp_length), 2)
            out[name] = str(completeness)

    # Output
    with open(f'{str(args.sample)}_amplicon_completeness.tsv', 'w') as f:
        header = 'sample\t{0}'.format('\t'.join(out.keys()))
        line = '{0}\t{1}'.format(str(args.sample),"\t".join(out.values()))
        f.write(header)
        f.write('\n')
        f.write(line)
        f.write('\n')

if __name__ == '__main__':
    main()
