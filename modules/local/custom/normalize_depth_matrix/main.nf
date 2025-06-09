process NORMALIZE_DEPTH_MATRIX {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pandas:1.5.2"

    input:
    path(depth_tsvs)

    output:
    path "sample_positional_normalized_depth.csv", emit: full_csv
    path "sample_positional_normalized_depth_stats.csv", emit: stats_csv

    script:
    """
    #!/usr/bin/env python3
    '''
    Simple script to normalize and merge together all positional depth
    '''
    import pandas as pd
    from pathlib import Path

    def col_max(df: pd.DataFrame, col: str) -> int:
        '''Find max depth of the given col in given df'''
        return max(df[col])

    # Summarize
    inpath = Path('./')
    outl = []
    first = True

    # The names are very much always the same so can use the '_depth' to get the name
    for bed in inpath.glob('*.tsv'):
        name = str(bed).split('_depth')[0]

        df = pd.read_csv(bed, sep='\t', names=['chrom', 'position', 'depth'])
        # Don't summarize empty files
        if df.empty:
            continue
        max_depth = col_max(df, 'depth')

        # Don't break on max depth of 0, each position is still 0
        #  the normalization is position/max and each position is 0
        if max_depth == 0:
            max_depth = 1

        df[name] = ((df['depth'] / max_depth) * 100).round(2)

        if first:
            outl.append(df[['position', name]])
            first = False
        else:
            outl.append(df[[name]])

    concat_df = pd.concat(outl, axis=1)
    concat_df = concat_df.set_index('position')
    out_df = concat_df.copy()
    out_df_t = out_df.transpose()
    out_df_t.to_csv('sample_positional_normalized_depth.csv')

    concat_df['mean'] = concat_df.mean(axis=1)
    concat_df['q1'] = concat_df.quantile(0.25, axis=1)
    concat_df['q3'] = concat_df.quantile(0.75, axis=1)
    concat_df['IQR'] = concat_df['q3'] - concat_df['q1']
    concat_df = concat_df[['mean', 'q1', 'q3', 'IQR']]
    concat_df.to_csv('sample_positional_normalized_depth_stats.csv')
    """

    stub:
    """
    touch sample_positional_normalized_depth.csv
    touch sample_positional_normalized_depth_stats.csv
    """
}
