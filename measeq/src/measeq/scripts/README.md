# Files
The files here are used within the `measeq` script to adjust how viral recon is run and carry out some of the pipeline's post-processing functions.

## [`Nextflow Measles Configuration File`](measles_parameters.config)
- The `measles_parameters.config` file adjusts some small aspects of the viralrecon pipeline to better handle the incoming measles Illumina sequencing data including:
    - Tool version updates
    - Tool parameter adjustments

## [`Rmarkdown Report File`](MeaSeq_Report.Rmd)
- The `MeaSeq_Report.Rmd` creates the MeaSeq Report output in final results directory using all of the data generated along with the additional Rmd sub-files.

## [`INDEL Check`](identify_indels.R)
- This R script is used to identify positions with the sample's genome where base insertions or deletions occured.

## [`Sample Variation`](calc_bam_variation.py)
- This Python script calculates the base variation by parsing the sample's BAM file reporting positions with base variation above specified thresholds (read count: 10, non-reference base perecntage: 15%).

## [`Logo: Public Health Agency of Canada`](phac.svg)
- The logo of the Public Health Agency of Canada is used within the [`Rmarkdown Report File`](Run_Report.Rmd) for branding purposes.