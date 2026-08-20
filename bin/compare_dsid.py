#!/usr/bin/env python3
'''Compare all sample N450 sequences to reference database and assign new ones a hash dsid'''
import argparse
import csv
import gzip
import hashlib
from Bio import SeqIO
from pathlib import Path

AMBIGUOUS_NT = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V']


def init_parser() -> argparse.ArgumentParser:
    """
    Specify command line arguments
    Returns command line parser with inputs
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f',
        '--fasta',
        required=True,
        type=Path,
        help='Path to input multifasta file to check dsids for'
    )
    parser.add_argument(
        '-d',
        '--dsid_fasta',
        required=False,
        type=Path,
        help='Path to input DSID Fasta database'
    )
    parser.add_argument(
        '-o',
        '--outfile',
        required=False,
        default='dsids.tsv',
        type=str,
        help='Output name for the final dsid results (Default: dsids.tsv)'
    )
    parser.add_argument(
        '-w',
        '--write_novel',
        action='store_true',
        help='Write novel sequences to output novel_dsids.tsv file'
    )

    return parser


def load_dsids(dsid_fasta: Path) -> dict:
    '''
    Purpose:
    --------
    Load DSIds into a dict to compare to

    Parameters:
    -----------
    dsid_fasta - Path
        Path to dsid fasta file

    Returns:
    --------
    Dict structured as {sequence: disd}
    '''
    which_open = gzip.open if dsid_fasta.suffix == '.gz' else open
    with which_open(dsid_fasta, 'rt') as f:
        dsid_seqs = {str(record.seq).upper(): record.id for record in SeqIO.parse(f, "fasta")}
    return dsid_seqs


def hash_seq(seq: str, length: int = 7) -> str:
    '''
    Purpose:
    --------
    Hash newly found N450 sequences so that they match between runs

    Parameters:
    -----------
    seq - str
        N450 sequence string
    length - int
        How long of the hash to return

    Returns:
    --------
    String of input length of the hash
    '''
    return hashlib.md5(seq.encode()).hexdigest()[:length]


def label_seqs(input_fasta: Path, dsid_seqs: dict, out: str, dsid_fasta: Path | None = None) -> dict:
    '''
    Purpose:
    --------
    Match and write detected IDs to TSV outfile along with tracking new IDs

    Parameters:
    -----------
    input_fasta - Path
        Path to input multifasta file
    dsid_seqs - dict
        Dict to match sequences to structured as {sequence: disd}
    out - str
        Output name for the files
    dsid_fasta - Path
        Path to dsid database fasta to extract file name

    Returns:
    --------
    Dict structured as {sequence: disd} for new sequences
    '''
    # For if novel seqs found
    novel_seqs = {}

    # Get the name of the dsid database fasta file without extensions
    if dsid_fasta:
        dsid_fasta_name = dsid_fasta.name
        for ext in ['.gz', '.fasta']:
            if dsid_fasta_name.endswith(ext):
                dsid_fasta_name = dsid_fasta_name.removesuffix(ext)
    else:
        dsid_fasta_name = 'None'

    # Match
    with open(out, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(['sample', 'matched_dsid', 'completeness', 'dsid_file_used'])

        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(record.seq).upper()
            sample = record.id.split('-N450')[0]

            # Check for N's and IUPACs before full matches
            completeness = 100.00
            if len(seq) == 0:
                dsid = 'No Data'
            elif 'N' in seq:
                n_pct = round((seq.count('N')/len(seq)) * 100, 2)
                completeness = 100 - n_pct
                dsid = 'Incomplete'
            elif any(nt in seq for nt in AMBIGUOUS_NT):
                dsid = 'Ambiguous Base'
            elif seq in dsid_seqs:
                dsid = dsid_seqs[seq]
            else:
                if seq not in novel_seqs:
                    seq_hash = hash_seq(seq)
                    label = f'Novel-{seq_hash}'
                    novel_seqs[seq] = label
                dsid = novel_seqs[seq]

            writer.writerow([sample, dsid, completeness, dsid_fasta_name])

    return novel_seqs


def main():
    '''Main entry point'''
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Create DSID Dict if we have data
    dsid_seqs = {}
    if args.dsid_fasta:
        dsid_seqs = load_dsids(args.dsid_fasta)

    # Label and output
    novel_seqs = label_seqs(args.fasta, dsid_seqs, args.outfile, args.dsid_fasta)

    if novel_seqs and args.write_novel:
        with open('novel_dsids.tsv', 'w', newline='') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            writer.writerow(['novel_id', 'sequence'])
            for seq, novel_id in novel_seqs.items():
                writer.writerow([novel_id, seq])

if __name__ == '__main__':
    main()
