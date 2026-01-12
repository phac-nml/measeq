# phac-nml/MeaSeq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.0.0] - 2026-01-12

Initial full pipeline release that includes equivalent Illumina and Nanopore workflows with full genome consensus sequence generation, N450 reporting, DSId hashing and assigning, and a final QC report.

This full release adds in support for all 24 genotypes when running with the genotype predictions provided the user sets up all the appropriate files. Currently: A, B3, and D8 are suppored in this repo although its recommended that users determine their own reference and primer files as they may not match the defaults.

### `Added`

- Support for all 24 genotypes and their primer files [PR #26](https://github.com/phac-nml/measeq/pull/26/files)

  - Recommended to set these with a `-params-file` if setting up multiple to make it easier to rerun the pipeline

- Reenabled support for contact information in [PR #26](https://github.com/phac-nml/measeq/pull/26/files) to be added to the final report using any combination of:
  - `--contact_name`: Name(s) to put on the contact page
  - `--contact_phone`: Phone number
  - `--contact_email`: Email
  - `--contact_website`: Website URL

### `Adjusted`

- Handling of intermediate files to allow the use of a full measles genome and/or genomes with the 5' and 3' UTR cut as reference

- Splitting of amplicon data to be per genome in case of the use of different amplicons across reference files

- Reorganization of the final report's mean genomic depth to be per genome and on its own tab in case of multiple references

## [v0.5.0] - 2025-12-08

The pipeline has been reorganized to run each sample with it's own reference to allow for the prediction of each sample's genotype and its mapping to its appropriate reference.

The default running inputs have been simplified to not require the `--reference` parameter. The parameter is still retained and able to be used by the user. However, by default, the pipeline will now automatically do genotyping and reference assignment. This is done to improve the accuracy of the outputs and streamline the workflow when there is a mix of sample genotypes.

### `Added`

- `PREDICT_GENOTYPE` is a python script that uses `minimap2`, `samtools`, and the WHO's N450 measles database to predict the most likely genotype of each sample and automatically selects the appropriate reference for downstream mapping [PR #24](https://github.com/phac-nml/measeq/pull/24)

- `EXTRACT_GENOTYPE` is a quick script that retries each reference's genotype and adds it the the reference's meta map for later downstream processes [PR #24](https://github.com/phac-nml/measeq/pull/24)

- New paramter for read filtering for ONT data [PR #24](https://github.com/phac-nml/measeq/pull/24)
  - `--ont_min_read_length`

### `Adjusted`

- Major changes to the pipeline's workflow to create a channel for every process input that requires matching sample and reference data to ensure each process uses the correct reference files that match the sample's genotype[PR #24](https://github.com/phac-nml/measeq/pull/24)

- Various minor tweaks have been made to input/output handling within processes across the pipeline to support the per-sample worklfow [PR #24](https://github.com/phac-nml/measeq/pull/24)

- The report generation process was adjusted to allow for a higher memory threshold to allow more samples to be included within one run [PR #24](https://github.com/phac-nml/measeq/pull/24)

- Tests have been updated and expanded to support the more complicated channel structure and better assess future updates [PR #24](https://github.com/phac-nml/measeq/pull/24)
  - CI workflow file reformat done to separate out the full workflow test done with the test profiles to be under the `ci_workflow.yml` file [PR #25](https://github.com/phac-nml/measeq/pull/25)
    - This is to prevent running the nf-test tests twice as much as necessary

## [v0.4.3] - 2025-11-07

Adjusting alignment filtering for both nanopore and Illumina data to remove supplementary reads and secondary reads. The supplementary reads were rarely adding in artifacts to the final consensus sequence including SNPs and INDELs that only they contained. This lowers read counts slightly but provides more accurate consensus sequences based on testing.

Also added in updates to the Clair3 process and Nanopore workflow including bumping the Clair3 tool version and adjusting the VCF filter found in the custom [Artic VCF Filter](https://github.com/artic-network/fieldbioinformatics/blob/master/artic/vcf_filter.py) to add in their SNP allele frequency and INDEL frameshift cutoff along with an addition of a RefCall filter. Note that the RefCall filter doesn't affect the consensus result at all, just the report mutation table.

Overall, the update should fix some rare artifacts to improve the accuracy of the final results

### `Added`

- `SAMTOOLS_SORT` nf-core module to replace the `SAMTOOLS_INDEX` module in Illumina pipeline [PR #21](https://github.com/phac-nml/measeq/pull/21)
  - As `sort` was removed from the process `BWAMEM2_MEM` itself to use Samtools view and filter specific alignments
  - Output remains the same as `<sample>.sorted.bam` with the alignment filtering having been done

### `Adjusted`

- `BWAMEM2_MEM` process adjusted to run Samtools view with `-f 3 -F 2048` to just keep mapped and properly paired alignments along with removing supplementary [PR #21](https://github.com/phac-nml/measeq/pull/21)
- `MINIMAP2_ALIGN` process added in the Samtools view flag `-F 2308` to remove unmapped, not primary alignment, and supplementary alignments [PR #21](https://github.com/phac-nml/measeq/pull/21)
- `Clair3` bumped to version 1.2.0 along with adjusting the command arguments based on internal sample testing
- Tests adjusted to conform with the overall lower mapped reads with the change in both variant calling subworkflows [PR #23](https://github.com/phac-nml/measeq/pull/23)
  - There was a single variant change in nanopore data test
    - Lower frequency variant at 0.47 AF
  - A couple more Ns overall at the start and end of genome for both
  - Slightly lower depth
- `nanoq` removed the max read length from `modules.json`
- `cs_vcf_filter.py` adjusted to match the source more with adding in min frameshift quality (30 instead of 50) and min allele freq [PR #21](https://github.com/phac-nml/measeq/pull/21)
  - Added a RefCall filter so that RefCall were just removed as well

## [v0.4.2] - 2025-10-23

Exposing more parameters to allow users more options in adjusting Illumina data analyses with parameters

### `Added`

- New parameters to control Illumina variant calling [PR #18](https://github.com/phac-nml/measeq/pull/18)
  - `--ivar_trim_min_read_length`
  - `--min_alt_fraction_freeabyes`
  - `--min_variant_qual_freebayes`
- Minimum and Maximum thresholds for parameters in the `nextflow_schema` [PR #18](https://github.com/phac-nml/measeq/pull/18)
- Example DSID fasta file [PR #19](https://github.com/phac-nml/measeq/pull/19)

### `Fixes`

- Fixed the final report provenance typo [PR #18](https://github.com/phac-nml/measeq/pull/18), [Issue #17](https://github.com/phac-nml/measeq/issues/17)

## [v0.4.1] - 2025-09-09

Small addition of Picard MarkDuplicates workflow along with some new tests

### `Added`

- nf-core Picard MarkDuplicates workflow as an optional parameter/workflow to use for Illumina data [PR #15](https://github.com/phac-nml/measeq/pull/15)
  - Along with this, added the bam stats samtools workflow to run even when the picard workflow isn't to keep outputs the same

## [v0.4.0] - 2025-09-03

### `Added`

- 3 columns to final Excel and CSV file [PR #12](https://github.com/phac-nml/measeq/pull/12)
  - N450_completeness
  - N450_mean_depth
  - N450_status
- Fail tracking for all samples and runs where every sample fails [PR #12](https://github.com/phac-nml/measeq/pull/12)

### `Adjusted`

- Adjusted the DSId final report page to have a list of all individual sample calls and the summary data shown [PR #12](https://github.com/phac-nml/measeq/pull/12)
- Added the N450 status to the initial report summary page [PR #12](https://github.com/phac-nml/measeq/pull/12)
- Always run the DSId check even with no database fasta file [PR #12](https://github.com/phac-nml/measeq/pull/12)
  - Everything will be labeled as `novel-hash` but it will group them up
- Updated CI tests [PR #12](https://github.com/phac-nml/measeq/pull/12)

## [v0.3.2] - 2025-08-01

### `Added`

- DSId tab to final HTML report for DSId summary information if a DSId file was given as input [PR #10](https://github.com/phac-nml/measeq/pull/10)
- New `overall.xlsx` final QC file based on adding in a few more columns [PR #10](https://github.com/phac-nml/measeq/pull/10)
  - `genome_fasta` and `N450_fasta` that contain fasta formatted sequence data
    - This broke the CSV file version parsing so that remains unchanged to get metadata to IRIDANext

### `Adjusted`

- Fixed a bug with MeaSeq Report summary table not correctly linking to certain samples [PR #10](https://github.com/phac-nml/measeq/pull/10)
- Removed `versions.yml` from being created during final report as versions already were reported [PR #10](https://github.com/phac-nml/measeq/pull/10)
- Adjusted code for final report row formatting to match throughout Rmd files [PR #10](https://github.com/phac-nml/measeq/pull/10)
- DSId assignment changes [PR #10](https://github.com/phac-nml/measeq/pull/10)
  - Hash adjusted to 7 characters
  - Semi-Complete removed to just be Incomplete
- New container for the MAKE_FINAL_QC_CSV step to remove artic container [PR #10](https://github.com/phac-nml/measeq/pull/10)
  - Adds in openpyxl

## [v0.3.1] - 2025-07-29

### `Adjusted`

- Fixed a bug where the irida json file wasn't being populated with metadata [PR #8](https://github.com/phac-nml/measeq/pull/8)
- Adjusted Illumina nf-test to use the `sample_name` field of the samplesheet to test that it works [PR #8](https://github.com/phac-nml/measeq/pull/8)

### `Removed`

- Normalized median read depth plot for now as it overlaps the full depth one too much [PR #8](https://github.com/phac-nml/measeq/pull/8)

## [v0.3.0] - 2025-07-25

### `Added`

- CI tests and Linting [PR #6](https://github.com/phac-nml/measeq/pull/6)
- Providence to final `Measeq_Report.html` file [PR #5](https://github.com/phac-nml/measeq/pull/5)
- Approximate genome position annotation to final report figures [PR #5](https://github.com/phac-nml/measeq/pull/5)
- `nf-validation` to work with IRIDA Next [PR #6](https://github.com/phac-nml/measeq/pull/6)
- `nf-iridanext` to work with IRIDA Next [PR #6](https://github.com/phac-nml/measeq/pull/6)
- `nf-prov` to allow some more providence options [PR #5](https://github.com/phac-nml/measeq/pull/5)
- Hashes to novel DSId calls so that in the same run new DSIds will match [PR #5](https://github.com/phac-nml/measeq/pull/6)
  - And they will match in other locations as well
- `min_indel_threshold` parameter to set the minimum indel threshold required to call an indel [PR #6](https://github.com/phac-nml/measeq/pull/6)

### `Adjusted`

- Output directories for Nanopore process in the `modules.config` file [PR #6](https://github.com/phac-nml/measeq/pull/6)
- Updated `artic` version from `1.6.2` --> `1.7.4` for nanopore pipeline [PR #6](https://github.com/phac-nml/measeq/pull/6)
- Ambiguous position handling for Illumina data [PR #5](https://github.com/phac-nml/measeq/pull/5)
  - Specifically for rare postions where there was a low-supported INDEL along with a SNP
- Negative control default string to add in 'en' [PR #5](https://github.com/phac-nml/measeq/pull/5)
- Samtools depth `meta1` to `meta` as it was breaking the IRIDA-Next plugin [PR #6](https://github.com/phac-nml/measeq/pull/6)
- Minimum nextflow version required to `24.10.0` [PR #5](https://github.com/phac-nml/measeq/pull/5)
- Samplesheet input added a `sample_name` column to work with IRIDA Next [PR #6](https://github.com/phac-nml/measeq/pull/6)

### `Removed`

- nf-schema and associated workflows [PR #6](https://github.com/phac-nml/measeq/pull/6)

## [v0.2.1-dev] - 2025-06-12

### `Added`

- New plots to the summary table for the final report
  - Summary sequencing depth plots
  - No longer stand-alone as it can get a bit large, a smaller report will be made later
- Allow primer schemes to have different direction extensions in the names
  - `_L`, `_R`, `_FORWARD`, `_REVERSE`, `_F`
  - Added in proper errors for if the primer file was not formatted correctly
- Readme updates for primers

### `Bugfixes`

- `genome_completeness` fixed to `genome_completeness_percent` in the negative control checking

### `Removed`

- N plots of the final report to help lower size and speed up opening

## [v0.2.0-dev] - 2025-06-06

### `Added`

- All of the pipeline has been rewritten in Nextflow
- Illumina paired-end sequencing workflow added
  - Freebayes for variant calling over ivar variants/consensus previously
- Nanopore (initial) workflow added
  - clair3
- DSId assignment added when using `--dsid fasta` parameter
  - Based on full sequence match
- Summary outputs added
  - Amplicon summary report
  - The current Rmarkdown report needs to be fixed for the new outputs

### `Deprecated`

- Current MeaSeq script that utilized viralrecon depreciated to make the whole pipeline nextflow

## v0.1.0 - 2025-04-08

### `Added`

- MeaSeq pipeline created and initial code added

[v1.0.0]: https://github.com/phac-nml/measeq/releases/tag/1.0.0
[v0.5.0]: https://github.com/phac-nml/measeq/releases/tag/0.5.0
[v0.4.3]: https://github.com/phac-nml/measeq/releases/tag/0.4.3
[v0.4.2]: https://github.com/phac-nml/measeq/releases/tag/0.4.2
[v0.4.1]: https://github.com/phac-nml/measeq/releases/tag/0.4.1
[v0.4.0]: https://github.com/phac-nml/measeq/releases/tag/0.4.0
[v0.3.2]: https://github.com/phac-nml/measeq/releases/tag/0.3.2
[v0.3.1]: https://github.com/phac-nml/measeq/releases/tag/0.3.1
[v0.3.0]: https://github.com/phac-nml/measeq/releases/tag/0.3.0
[v0.2.1-dev]: https://github.com/phac-nml/measeq/releases/tag/0.2.1-dev
[v0.2.0-dev]: https://github.com/phac-nml/measeq/releases/tag/0.2.0-dev
