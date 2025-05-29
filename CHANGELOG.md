# phac-nml/MeaSeq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.2.0 - [2025-05-22]

### `Added`

- All of the pipeline has been rewritten in Nextflow
- Illumina paired-end sequencing workflow added
    - Freebayes for variant calling over ivar variants/consensus previously
- Nanopore (initial) workflow added
    - clair3
- DSID assignment added when using `--dsid fasta` parameter
    - Based on full sequence match
- Summary outputs added
    - Amplicon summary report
    - The current Rmarkdown report needs to be fixed for the new outputs

### `Deprecated`

- Current MeaSeq script that utilized viralrecon depreciated to make the whole pipeline nextflow

## v0.1.0 - [2025-04-08]

### `Added`

- MeaSeq pipeline created and initial code added
