#!/usr/bin/env python3
'''Create Summary QC CSV file based on pipeline outputs'''

import argparse
import pandas as pd

def init_parser() -> argparse.ArgumentParser:
    """
    Specify command line arguments
    Returns command line parser with inputs
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c',
        '--csv',
        required=True,
        type=str,
        help='Concatenated sample QC files'
    )
    parser.add_argument(
        '-m',
        '--metadata',
        required=False,
        type=str,
        help='Metadata TSV file to add to failed samples'
    )
    parser.add_argument(
        '--version',
        required=True,
        type=str,
        help='Pipeline version to include in output'
    )
    parser.add_argument(
        '-F',
        '--fill_str',
        required=False,
        default='NA',
        type=str,
        help='What to fill in columns with missing data. Default is "NA"'
    )
    parser.add_argument(
        '-T',
        '--threshold',
        required=False,
        default=10.0,
        type=float,
        help='Proportion genome completeness threshold to fail negative control samples by. Default: 10.0'
    )
    parser.add_argument(
        '--neg_ctrl_substrings',
        required=False,
        default='neg,ntc,blank',
        type=str,
        help='Comma separated substrings to be used to find negative control samples. Default: "neg,ntc,blank"'
    )
    return parser

def validate_df_columns(df: pd.DataFrame, needed_columns: list) -> None:
    """
    Purpose
    -------
    Check that input CSV contains the correct columns needed. Exits program if not

    Parameters
    ----------
    df: pd.DataFrame
        Pandas dataframe made from the input CSV file
    """
    columns = list(df.columns)
    if any(x not in columns for x in needed_columns):
        raise ValueError('Missing {} column(s) needed for validation'.format([x for x in needed_columns if x not in columns]))

def assess_control(row: pd.Series, threshold: float) -> str:
    """Assess control values to pass or fail them"""
    if row['genome_completeness'] >= threshold:
        return f'Warning - Above {threshold}% genome completeness contamination threshold'
    return 'PASS'

def main() -> None:
    '''Run the program'''
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()
    neg_ctrl_substrings = args.neg_ctrl_substrings.split(',')

    # Fill these columns with 0 if they have no data
    numeric_columns = [
        'num_input_reads',
        'num_aligned_reads',
        'num_consensus_n',
        'genome_completeness_percent',
        'mean_sequencing_depth',
        'median_sequencing_depth',
        'total_variants',
        'num_snps',
        'num_deletions',
        'num_deletion_sites',
        'num_insertions',
        'num_insertion_sites',
        'genome_length'
    ]

    # Do stuff
    df = pd.read_csv(args.csv)
    validate_df_columns(df, ['sample', 'num_aligned_reads', 'genome_completeness_percent', 'mean_sequencing_depth', 'median_sequencing_depth', 'divisible_by_6', 'qc_status'])

    # Add in metadata if there is any
    if args.metadata:
        metadata_df = pd.read_csv(args.metadata, sep='\t')
        validate_df_columns(df, ['sample'])
        df = df.merge(metadata_df, on='sample', how='left')

    # Check negative controls, can add more here
    neg_df = df[df['sample'].str.contains('|'.join(neg_ctrl_substrings), na=False, case=False)]
    if not neg_df.empty:
        # Check neg columns
        run_control_status = 'PASS'
        run_control_info = ''

        # Apply checks on controls and update normal df qc_status column
        df['neg_control_info'] = neg_df.apply(assess_control, threshold=args.threshold, axis=1)
        df.loc[df['neg_control_info'] == 'PASS', 'qc_status'] = 'PASS'
        df.loc[df['neg_control_info'].str.contains('Warning', na=False), 'qc_status'] = 'CONTROL_WARN'

        # Have to adjust this later
        if any(neg_df['genome_completeness'] >= args.threshold):
            failing_samples = neg_df[neg_df['genome_completeness'] >= args.threshold]['sample'].to_list()
            run_control_status = 'WARN'
            run_control_info = f'Samples: {";".join(failing_samples)} are above {args.threshold} contamination threshold'
    else:
        # Add neg control columns as not available
        run_control_status = 'WARN'
        run_control_info = 'No negative controls found in run'

    # Adding final columns and output
    df['run_status'] = run_control_status
    df['run_summary'] = run_control_info
    df = df.fillna(args.fill_str)
    df['pipeline_name'] = 'measeq'
    df['pipeline_version'] = args.version
    df.sort_values(by='sample', inplace=True)
    df.to_csv('overall.qc.csv', index=False)

if __name__ == '__main__':
    main()
