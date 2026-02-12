process AMPLICON_DEPTH_HEATMAP {
    label 'process_single'
    tag "$ref_id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:1.5.2'
        : 'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(ref_id), path(amplicon_beds)

    output:
    tuple val(ref_id), path ("${ref_id}_amplicon_depth_heatmap_mqc.tsv"), emit: heatmap_tsv
    path("${ref_id}_amplicon_depth_full.tsv"), emit: full_tsv

    script:
    """
    #!/usr/bin/env python3
    '''
    Simple script to merge together, transpose, and log10 all of the amplicon coverage files
    '''
    import glob
    import pandas as pd
    import numpy as np
    from pathlib import Path

    # Find files and take only the two needed columns
    tsv_files = glob.glob('*_amplicon_stats.bed')
    df_list = []
    for f in tsv_files:
        f = Path(f)
        name = f.name.split('_amplicon_stats.bed')[0]
        df = pd.read_csv(f, sep='\t', header=None)
        df = df[[3, 6]]
        df.columns = ['amplicon_id', name]
        df_list.append(df)

    # Transpose and Output
    df = pd.DataFrame().join([d.set_index('amplicon_id') for d in df_list], how='outer').dropna().reset_index()
    df.set_index('amplicon_id', inplace=True, drop=True)
    df = df.transpose()
    df.index.name = 'sample'
    df.sort_index(inplace=True)
    df.to_csv('${ref_id}_amplicon_depth_full.tsv', index=True, sep="\t")

    # Log10
    df1 = pd.DataFrame(np.ma.log10(df.values).filled(0), index=df.index, columns=df.columns)
    df1.to_csv('${ref_id}_amplicon_depth_heatmap_mqc.tsv', index=True, sep="\t")
    """

    stub:
    """
    touch ${ref_id}_amplicon_depth_full.tsv
    touch ${ref_id}_amplicon_depth_heatmap_mqc.tsv
    """
}
