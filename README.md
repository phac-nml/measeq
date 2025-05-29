# MeaSeq: Measles Sequence Analysis Automation

- [Updates](#updates)
  - [2025-05-29](#2025-05-29)
- [Introduction](#introduction)
- [Installation](#installation)
- [Resource Requirements](#resources-requirements)
- [Usage](#usage)
  - [Illumina](#illumina)
  - [Nanopore](#nanopore)
  - [Amplicon and Primer Files](#amplicon-and-primer-files)
  - [More Run Options](#more-run-options)
  - [Testing](#testing)
- [Outputs](#outputs)
- [Steps](#steps)
  - [Illumina Steps](#illumina-steps)
  - [Nanopore Steps](#nanopore-steps)
- [Troubleshooting](#troubleshooting)
- [Credits](#credits)
- [Citations](#citations)
- [Contributing](#legal)
- [Legal](#legal)

## Updates

### *2025-05-29*

- Switched to running all steps with nextflow for the following reasons:
  - Allow more control over all of the steps
  - Easier to install/run along with having more dependency management options (IE not required to use `conda`)
  - Eventual implementation to IRIDA-Next
- Focus is **currently on Illumina data** although the nanopore side *should* still work

## Introduction

**MeaSeq** is a measles virus (MeV) specific pipeline established for use in surveillance and outbreak analysis. This pipeline utilizes a reference-based read mapping approach for Whole Genome or Amplicon sequencing data from both the Illumina and Nanopore platforms to output MeV consensus sequences, variant data, sequencing qualtiy information, and custom summary reports.

![MeaSeq Workflow Diagram](todo)

This project aims to implement an open-source, easy to run, MeV Whole Genome Sequence analysis pipeline that works on both Illumina and Nanopore data. The end goal of this project is to deploy a standardized pipeline focused on final reporting metrics and plots for rapid detection and response to MeV outbreaks in Canada and abroad.

## Installation

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Installation requires both [nextflow](https://www.nextflow.io/) at a minimum version of `24.04.2` and a dependency management system to run.

Steps:

1. Download and install nextflow

   1. Download and install with [conda](https://docs.conda.io/en/latest/miniconda.html)
      - Conda command: `conda create -n nextflow -c conda-forge -c bioconda nextflow`
   2. Install with the instructions at https://www.nextflow.io/

2. Determine which dependency management system works best for you

   - _Note_: Currently the plotting process is using a custom docker container but it should work for both docker and singularity

3. Run the pipeline with one of the following profiles to handle dependencies (or use your [own profile](https://nf-co.re/docs/usage/getting_started/configuration) if you have one for your institution!:
   - `conda`
   - `mamba`
   - `singularity`
   - `docker`

## Resources Requirements

By default, the `bwamem2` step has a minimum resource usage allocation set to `12 cpus` and `72GB memory` using the nf-core `process_high` label.

This can be adjusted (along with the other labels) by creating and passing a [custom configuration file](https://nf-co.re/docs/usage/getting_started/configuration) with `-c <config>`. More info can be found in the [usage doc](./docs/usage.md)

The pipeline has also been test using as low as `2 cpus` and `16GB memory`

## Usage

### Illumina

First, prepare a samplesheet with your input data that looks as follows for Illumina paired-end data:

**samplesheet.csv**:

```csv
sample,fastq_1,fastq_2
MeVSample01,/PATH/TO/inputread1_S1_L002_R1_001.fastq.gz,/PATH/TO/inputread1_S1_L002_R2_001.fastq.gz
PosCtrl01,/PATH/TO/inputread2_S1_L003_R1_001.fastq.gz,/PATH/TO/inputread2_S1_L003_R2_001.fastq.gz
Sample3,/PATH/TO/inputread3_S1_L004_R1_001.fastq.gz,/PATH/TO/inputread3_S1_L004_R2_001.fastq.gz
```

Each row represents a sample and its associated paired-end Illumina read data.

You can then run the pipeline using:

```bash
nextflow run phac-nml/measeq \
    --input <SAMPLESHEET> \
    --outdir <OUTDIR> \
    --reference <REFERENCE FASTA> \
    --platform <illumina||nanopore> \
    -profile <docker/singularity/.../institute>
```

### Nanopore

And as follows for nanopore data:

**samplesheet.csv**

```csv
sample,fastq_1,fastq_2
MeVSample01,/PATH/TO/inputread1.fastq.gz,
PosCtrl01,/PATH/TO/inputread2.fastq.gz,
Sample3,/PATH/TO/inputread3.fastq.gz,
```

Each row represents a sample and its single-end nanopore data.

You can then run the pipeline using:

```bash
nextflow run phac-nml/measeq \
    --input <SAMPLESHEET> \
    --outdir <OUTDIR> \
    --reference <REFERENCE FASTA> \
    --platform <illumina||nanopore> \
    --clair3_model <MODEL> \
    -profile <docker/singularity/institute/etc>
```

### Amplicon and Primer Files

Both Illumina and Nanopore support running amplicon data using a primer scheme file. To run amplicon data all you need is a primer bed file where the primers have been mapped to the location in the reference genome used. The parameter being `--primer_bed <PRIMER_BED>`. An example primer bed file looks as such:

**primer.bed**

```
<CHROM>         <START> <END>   <PRIMER_NAME>   <POOL>  <DIRECTION>
MH356245.1      1       25      MSV_1_LEFT      1       +
MH356245.1      400     425     MSV_2_LEFT      2       +
MH356245.1      500     525     MSV_1_RIGHT     1       -
MH356245.1      900     925     MSV_2_RIGHT     2       -
```

_Note_: The first line is just to display what each line expects.

### More Run Options

For more detailed running options please refer to [the usage docs](docs/usage.md).

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

### Testing

To test the `MeaSeq` pipeline, and that everything works on your system, a small set of D8 genotype samples have been included from [SRA BioProject PRJNA480551](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA480551) in the [`test_data/fastqs`](test_dataset/fastqs/) directory.

To run the pipeline on these samples run the following command:

```bash
nextflow run phac-nml/measeq -profile test,<docker/singularity/institute/etc>
```

## Outputs

The main outputs of the pipeline are the `consensus sequences` (N450 and Full), the `overall.qc.csv` summary table, and the `MeaSeq_Report.html`. The final MeaSeq report gives a summary of the run including sample quality metrics, plots, and any additional information. Detailed pipeline outputs are described [within the output docs](docs/output.md)

## Steps

### Illumina Steps

1. Generate Reference and Primer Intermediates
2. FastQC
3. Illumina Consensus Workflow
    1. FastP
    2. BWAMem2
    3. Ivar Trim (Amplicon input only)
    4. Freebayes
    5. Process Freebayes VCF
    6. Make Depth Mask
    7. Bcftools Consensus (Ambiguous and Consensus variants)
4. Nextclade (N450 and Custom datasets, N450 fasta output)
5. Samtools depth
6. Compare DSID (Optional with `--dsid_fasta` parameter)
7. Make sample QC
8. Amplicon Summary Workflow (Amp only data)
    1. Bedtools Coverage
    2. Summarize Amplicon Depth
    3. Summarize Amplicon Completeness
    4. MultiQC Amplicon Report
9. Report Workflow
    1. Samtools mpileup
    2. Pysamstats
    3. Rmarkdown

### Nanopore Steps

To come

## Troubleshooting

For troubleshooting, please open an issue or consult [the usage docs](docs/usage.md) to see if they have the information you require.

## Credits

Written by McMaster University Co-op student - `Ahmed Abdalla`

For questions please contact either:

- Darian Hole (`darian.hole@phac-aspc.gc.ca`)
- Molly Pratt (`molly.pratt@phac-aspc.gc.ca`)

## Citations

> A citation for this pipeline will be available soon.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
> In addition, references of tools and data used in this pipeline are as follows:

Detailed citations for utilized tools are found in [CITATIONS.md](./CITATIONS.md)

## Contributing

Contributions are welcome through creating PRs or Issues

## Legal

Copyright 2025 Government of Canada

Licensed under the MIT License (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

https://opensource.org/license/mit/

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
