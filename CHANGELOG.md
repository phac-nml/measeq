# phac-nml/MeaSeq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.3.2] - 2026-09-01

Tiny adjustments including fixing the biopython conda environment, creating a channel for the `--ivar_primer_pairs` parameter and fixing a small reporting bug for low depth deletions that ultimately get ignored

### `Adjusted`

- Creating a chanel for `--ivar_primer_pairs` so they can be found in containers [PR #51](https://github.com/phac-nml/measeq/pull/51)
- Conda channel for biopython corrected to `conda-forge` [PR #53](https://github.com/phac-nml/measeq/pull/53)
- Rare issue where variant counts contained masked sites when the depth was below the depth threshold [PR #54](https://github.com/phac-nml/measeq/pull/54)

## [v1.3.1] - 2026-08-20

Allowing the `--dsid_fasta` file to be gzipped for input, adding in an internal specific default DSId file, and adding the dsid file name to the sample's metadata

### `Added`

- New `--default_gsp_dsid` which is just for internal usage but it sets a fall back default file for the DSId assignment [PR #47](https://github.com/phac-nml/measeq/pull/47)
- New metadata column `dsid_file_used` for better tracking of which DSId database version the sample was compared against [PR #49](https://github.com/phac-nml/measeq/pull/49)
- Final report DSId summary page now contains value boxes with the name of the DSId file the samples were compared against and the most prevalent DSId [PR #49](https://github.com/phac-nml/measeq/pull/49)

### `Adjusted`

- DSId fasta file can now be gzipped for input if wanted [PR #45](https://github.com/phac-nml/measeq/pull/45)
- Processes that used stdout as an output were modified to use environment variables to fix an issue where the ouput had a newline within IRIDA Next [PR #48](https://github.com/phac-nml/measeq/pull/48)
- While loops were removed from pipeline setup as they were not supported by nextflow [PR #48](https://github.com/phac-nml/measeq/pull/48)
- Fixed an issue where sample names didn't link correctly in the final report [PR #48](https://github.com/phac-nml/measeq/pull/48)

## [v1.3.0] - 2026-07-13

Code cleanup, exposed nanopore frameshift quality parameter, added excess ambiguity warning, artic version bump to `v1.10.3`, and clair3 version bump from `v1.2.0` to `v2.0.2`.

Important Note: The Clair3 `v1.2.0` models have to be converted to the new format. The update to `artic` will pull these new versions by default if you are providing a `--model <model name>` name on the command line. If you are using a local model you'll have to grab the new version. More info is available in [clair3's docs](https://github.com/HKU-BAL/clair3#pre-trained-models). No consensus variant output changes noted in the test datasets with this update.

### `Added`

- `--min_frameshift_qual_c3` arg to adjust the minimum cutoff for a frameshift variant to be called for nanopore data [PR #43](https://github.com/phac-nml/measeq/pull/43)

  - Default is the same at `30`

- Added an excess ambiguity warning [PR #43](https://github.com/phac-nml/measeq/pull/43)

  - `>=5` sites for Illumina by default
  - `>=10` sites for Nanopore by default
  - Its configurable in the `modules.json` file

- Added nextflow channel check to make sure that the input model folder has a `.pt` extension file required [PR #43](https://github.com/phac-nml/measeq/pull/43)

### `Adjusted`

- Renamed the VCF parsing scripts to make them easier to distinguish [PR #43](https://github.com/phac-nml/measeq/pull/43)

  - process_nanopore_vcf.py
  - process_illumina_vcf.py

- Illumina VCF parsing adjusted to better handle complex sites where an INDEL and a SNP are called with the SNP being the final call [PR #43](https://github.com/phac-nml/measeq/pull/43)

  - No consensus output differences
  - Adjusts the calculation for multi-SNP sites on what percentage each position is the reference or alt based on the CIGAR string of the variant
    - For IUPAC base selection
  - Adjusts the need to split ambiguous sites from consensus sites in VCF outputs for downstream processing
    - New filename for the illumina consensus variants file `.consensus.norm.vcf.gz` instead of `.processed.norm.vcf.gz`

- Split out the `artic` subcommands into their own folders to better match the nextflow best practices [PR #43](https://github.com/phac-nml/measeq/pull/43)
- Artic updated to `v1.10.3` in all spots that were `v1.8.5` [PR #43](https://github.com/phac-nml/measeq/pull/43)
- Clair3 updated to `v2.0.2` [PR #43](https://github.com/phac-nml/measeq/pull/43)

- Resource label adjustments [PR #44](https://github.com/phac-nml/measeq/pull/44):
  - `CALCULATE_BAM_VARIATION` moved to `process_medium` from `process_high_memory` as that isn't a standard label
  - Moving nf-core modules from `process_high` to `process_medium` as with the measles genome they don't require as many resources as we are giving
    - `BWAMEM2_MEM`
    - `BOWTIE2_BUILD`
    - `BOWTIE2_ALIGN`

## Removed

- Intermediate BCFTools consensus step in illumina workflow that applied ambiguous sites [PR #43](https://github.com/phac-nml/measeq/pull/43)
  - Ambiguous sites now have correct genotypes associated and can just run the command once

## [v1.2.0] - 2026-05-20

Update to enable users to map reads with Bowtie 2 instead of BWAMem as an optional parameter and use Artic primers as default. Add in parameter to better control ONT variant masking

### `Added`

- Illumina reads mapping with Bowtie 2 as an alterative (instead of BWAMem 2) in response to [issue #37](https://github.com/phac-nml/measeq/issues/37). [PR #39](https://github.com/phac-nml/measeq/pull/39)

  - `align_bowtie2` added as an optional parameter to map illumina data with Bowtie 2 instead of BWAMem2.

- [Artic primers](https://doi.org/10.1101/2024.12.20.629611) for MeV were added as a profile. [PR #39](https://github.com/phac-nml/measeq/pull/39)

  - This allows running the pipeline with the Artic primers mapped to the pipeline's preset references (D8, B3, and A genotypes) [issue #38](https://github.com/phac-nml/measeq/issues/38).
  - To run the pipeline with this profile, use `nextflow run phac-nml/measeq -profile artic_primers,<docker/singularity>` with the other normal parameters you would use.

- New `min_mask_freq_c3` parameter for nanopore data to help better control when sites are masked as Ns in more noisy regions/datasets [PR #42](https://github.com/phac-nml/measeq/pull/42)
  - Default is `0.30`
  - To help allow adjustments to the acceptable amount of site variation before an N is used to mask it
    - So now by default a call with `0.15` alt frequeny and a quality score of 3 would NOT be N masked when before it was

> [!NOTE]
> If you use `-profile artic_primers`, then there is no need to use `--amplicon` as it is automatically passed.

### `Fixes`

- Ambiguous regions that don't include the Ref base are properly tracked now in the report and include both alts [PR #41](https://github.com/phac-nml/measeq/pull/41)

## [v1.1.0] - 2026-03-27

Update focusing on fixing small inconsistencies between final reports and output consensus sequences along with a few small bugfixes and exposing some more parameter options available to the user

### `Fixes`

- Artic `align_trim` downgraded to `1.7.4` to fix issue where final position reads are being trimmed too strictly causing the last few bases to always be dropped [PR #33](https://github.com/phac-nml/measeq/pull/33)
- Depth calculations adjusted so that calculations between depth mask and samtools depth match in final report [PR #33](https://github.com/phac-nml/measeq/pull/33)
- Some nanopore variants were missed at the end of the genome [PR #33](https://github.com/phac-nml/measeq/pull/33)
  - Artic version in subcommands bumped to `v1.8.5` to solve this
  - Added a wrapper around `align_trim` to set that random seed at 42 [PR #35](https://github.com/phac-nml/measeq/pull/35)

### `Adjusted`

- Nanopore variant filter added a depth percentage check at a required minimum 5% of the total positional depth to make sure variants in overlapping regions are not extremely low in depth compared to the other amplicon [PR #33](https://github.com/phac-nml/measeq/pull/33)

  - If there is a massive discrepency and there is a true variant it should be in both so it will still be kept
  - If there is no variant in the other overlapping amplicon that would suggest an error so ignore
  - Parameter added called `min_site_threshold_c3` to handle this
  - LowQual sites ignored from masking

- Nanopore variant filter parameter `min_allele_freq_c3` added to allow more control to users [PR #33](https://github.com/phac-nml/measeq/pull/33)

- Multi-allelic illumina variants are labelled as such if they are lower depth and have been split up [PR #33](https://github.com/phac-nml/measeq/pull/33)

- New column for Nanopore data called `num_masked` which replaces `num_iupac` [PR #35](https://github.com/phac-nml/measeq/pull/35)
  - This was done we are not calling mixed sites as IUPAC bases in ONT data
  - It can be seen in the final CSV file, excel file, and final report

### `Added`

- Illumina primer trimming parameter additions [PR #33](https://github.com/phac-nml/measeq/pull/33)
  - `iVar trim` offset added as a parameter, set to 0 by default to help trim tricky primer combinations with parameter `ivar_offset`
  - `iVar trim` primers TSV file option added as an experimental hidden parameter called `ivar_primer_pairs` for if users who would like to use it

## [v1.0.1] - 2026-02-12

Adjusting containers to add in the `task.ext.override_configured_container_registry` to processes with containers not from quay. This allows more flexibility and the ability to use a custom or private container registry for all processes. Also meets our [IRIDA Next Criteria](https://github.com/phac-nml/pipeline-standards?tab=readme-ov-file#622-configuring-module-software-with-private-or-alternate-container-registries)

Updating the script that trims nanopore alignment ends to a newer version that takes into account a seed, ensuring the reproducibility of results across multiple runs

### `Adjusted`

- Processes where the container registry is not quay have `task.ext.override_configured_container_registry` added [PR #30](https://github.com/phac-nml/measeq/pull/30)
- Process `MAKE_CUSTOM_REPORT` adjusted from label `process_high_memory` to `process_medium` [PR #30](https://github.com/phac-nml/measeq/pull/30)
- Process `ARTIC_ALIGN_TRIM` adjusted to a newer version of ARTIC [PR#31](https://github.com/phac-nml/measeq/pull/31) and renamed `ALIGN_TRIM`

## [v1.0.0] - 2026-02-04

Initial full pipeline release that includes equivalent Illumina and Nanopore workflows with full genome consensus sequence generation, N450 reporting, DSId hashing and assigning, and a final QC report.

This full release adds in support for all 24 genotypes when running with the genotype predictions provided the user sets up all the appropriate files. Currently: A, B3, and D8 are suppored in this repo although its recommended that users determine their own reference and primer files as they may not match the defaults.

### `Added`

- Support for all 24 genotypes and their primer files [PR #26](https://github.com/phac-nml/measeq/pull/26/files)

  - Recommended to set these with a `-params-file` if setting up multiple to make it easier to rerun the pipeline

- New parameters to help allow more control to ONT data

  - `ont_min_base_qual`: The minimum base quality of ONT reads to keep
  - `ont_keep_incorrect_primers`: When running with amplicon data, keep reads that don't match up to their primer pair

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

[v1.2.0]: https://github.com/phac-nml/measeq/releases/tag/1.2.0
[v1.1.0]: https://github.com/phac-nml/measeq/releases/tag/1.1.0
[v1.0.1]: https://github.com/phac-nml/measeq/releases/tag/1.0.1
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
