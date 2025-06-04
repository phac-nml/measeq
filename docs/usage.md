# MeaSeq pipeline: Usage

## Introduction

This pipeline is intended to be run on measles virus (MeV) paired-end Illumina or single-end Nanopore sequencing data. It is written in nextflow to process MeV specific data with various outputs. This pipeline is intended for rapid deployment in outbreak situations in Canada and abroad.

## Index

- [Samplesheet Input](#samplesheet-input)
  - [Illumina samplesheet.csv](#illumina-samplesheetcsv)
  - [Nanopore samplesheet.csv](#nanopore-samplesheetcsv)
  - [Samplesheet Column Descriptions](#samplesheet-column-descriptions)
- [Running the Pipeline](#running-the-pipeline)
  - [Illumina Required](#illumina-required)
  - [Nanopore Required](#nanopore-required)
  - [Expanded Parameter Options](#expanded-parameter-options)
    - [Metadata TSV](#metadata-tsv)
    - [Clair3 Models](#clair3-models)
    - [DSId Matching](#dsid-matching)
    - [All Parameters Table](#all-parameters-table)
  - [Other Settings and Parameter Files](#other-settings-and-parameter-files)
  - [Updating the Pipeline](#updating-the-pipeline)
  - [Reproducibility](#reproducibility)
- [Core Nextflow Arguments](#core-nextflow-arguments)
  - [-profile](#-profile)
  - [-resume](#-resume)
  - [-c](#-c)
- [Custom Configuration](#custom-configuration)
  - [Resource Requests](#resource-requests)
  - [Custom Containers](#custom-configuration)
  - [Custom Tool Arguments](#custom-tool-arguments)
- [Running in the Background](#running-in-the-background)
- [Nextflow Memory Requirements](#nextflow-memory-requirements)

## Samplesheet Input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

### Illumina samplesheet.csv

```csv
sample,fastq_1,fastq_2
MeVSample01,/PATH/TO/inputread1_S1_L002_R1_001.fastq.gz,/PATH/TO/inputread1_S1_L002_R2_001.fastq.gz
PosCtrl01,/PATH/TO/inputread2_S1_L003_R1_001.fastq.gz,/PATH/TO/inputread2_S1_L003_R2_001.fastq.gz
Sample3,/PATH/TO/inputread3_S1_L004_R1_001.fastq.gz,/PATH/TO/inputread3_S1_L004_R2_001.fastq.gz
```

### Nanopore samplesheet.csv

```csv
sample,fastq_1,fastq_2
MeVSample01,/PATH/TO/inputread1.fastq.gz,
PosCtrl01,/PATH/TO/inputread2.fastq.gz,
Sample3,/PATH/TO/inputread3.fastq.gz,
```

Input Parameter:

```bash
--input '</PATH/TO/samplesheet.csv>'
```

### Samplesheet Column Descriptions

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

### Illumina Required

The typical base command for running the illumina pipeline path is as follows:

```bash
nextflow run phac-nml/measeq --input ./samplesheet.csv --outdir ./results  --reference ./REFERENCE.fa --platform illumina -profile docker
```

This will launch the pipeline with the `docker` configuration profile using `REFERENCE.fa` and for `illumina` input data. See [All Parameters](#all-parameters) for more information about all parameters available.

### Nanopore Required

The typical base command for running the nanopore pipeline path is as follows:

```bash
nextflow run phac-nml/measeq --input ./samplesheet.csv --outdir ./results  --reference ./REFERENCE.fa --platform nanopore --model 'r941_prom_sup_g5014' -profile docker
```

This will launch the pipeline with the `docker` configuration profile using `REFERENCE.fa` and for `nanopore` input data. We've also set the model for clair3 to be `r941_prom_sup_g5014`. See [All Parameters](#all-parameters) for more information about all parameters available.

### Expanded Parameter Options

Additional options to help run the pipeline to suit your needs

#### Metadata TSV

Metadata can be incorporated into the pipeline provided it is specified in TSV format with at minimum a `sample` column that matches to the samplesheet.csv input sample column

Example:

```tsv
sample	collection_date	etc.
MeV01	2024-05-10	...
MeV02	2024-06-12	...
MeV03	2024-09-05	...
```

An example file can be [found here](../assets/metadata.tsv)

#### Clair3 Models

The Nanopore pipeline utilizes [Clair3](https://github.com/HKU-BAL/Clair3) to call variants which requires a model that should be picked based off of the flowcell, pore, translocation speed, and basecalling model.

Some models are built into clair3 and some need to be downloaded. The [pre-trained clair3](https://github.com/HKU-BAL/Clair3?tab=readme-ov-file#pre-trained-models) models are able to be automatically downloaded when running the pipeline using [`artic get_models`](https://github.com/artic-network/fieldbioinformatics/blob/master/artic/get_models.py) and can be specified as a parameter with `--model <MODEL>`.

Additional or local models can also be used, you just have to provide a path to them and use the `--local_model <PATH>` parameter instead

#### DSId Matching

While 24 MeV genotypes were initially identified, only 2 have been detected since 2021: B3 and D8. Due to this, the Distinct Sequence Identifier (DSId) system was created to designate a unique 4-digit identifier based on the precise N450 sequence as a sub-genotype nomenclature. With this, the N450 sequence should be submitted to the [MeaNS2 database](https://who-gmrln.org/means2) to obtain one.

As this is a submission, a quick check is available to determine if you have any previously seen DSId's. To do this you have to setup a multifasta file with the DSId's as the header to match to the output N450 sequences and pass it in with `--dsid_fasta <MULTIFASTA>`. As an example:

```
> 1231
GTCAGTTCCACATTGGCATCT...
> 4412
GTCAGTTCCACATTGGCATCT...
> 5721
GTCAGTTCCACATTGGCATCT...
> etc.
GTCAGTTCCACATTGGCATCT...
```

#### All Parameters Table

A table containing all of the parameter descriptions. You can also do `nextflow run phac-nml/measeq --help` to get them on the command line

| Parameter                   | Description                                                                         | Required      | Type    | Default       | Notes                                            |
| --------------------------- | ----------------------------------------------------------------------------------- | ------------- | ------- | ------------- | ------------------------------------------------ |
| --input                     | Path to comma-separated file containing sample and read information                 | True          | Path    | null          |                                                  |
| --outdir                    | Name of output directory to store results                                           | True          | String  | null          |                                                  |
| --reference                 | Path to reference fasta file to map to                                              | True          | Path    | null          |                                                  |
| --platform                  | Sequencing platform used, either 'illumina or nanopore'                             | True          | Choice  | null          |                                                  |
| --model                     | Name of clair3 model to use                                                         | Nanopore data | String  | null          | Can use `--local_model` instead                  |
| --local_model               | Path to local clair3 model to use                                                   | Nanopore data | Path    | null          | Can use `--model` instead if wanted              |
| --primer_bed                | Path to bed file containing genomic primer locations                                | False         | Path    | null          | Use for amplicon data                            |
| --min_ambiguity_threshold   | Minimum threshold to call a position as an IUPAC                                    | False         | Float   | 0.30          | Illumina only                                    |
| --max_ambiguity_threshold   | Maximum threshold to call a position as an IUPAC                                    | False         | Float   | 0.75          | Illumina only                                    |
| --normalise_ont             | Normalise each amplicon barcode to set depth                                        | False         | Int     | 2000          | Nanopore only                                    |
| --min_variant_qual_c3       | Minimum variant quality to pass clair3 filters                                      | False         | Int     | 8             | Nanopore only                                    |
| --metadata                  | Path to metadata TSV file containing at minimum 'sample' column                     | False         | Path    | null          | See [Metadata TSV](#metadata-tsv)                |
| --dsid_fasta                | Path to DSID multi-fasta to match output consensus data to                          | False         | Path    | null          | See [DSId Matching](#dsid-matching)              |
| --min_depth                 | Minimum depth to call a base                                                        | False         | Int     | 10            |                                                  |
| --no_frameshifts            | Fail all indel variants not divisible by 3                                          | False         | Boolean | False         | Somewhat crude filter, only use if really needed |
| --neg_control_pct_threshold | Threshold of genome to be called in a negative control to fail it                   | False         | Int     | 10            |                                                  |
| --neg_ctrl_substrings       | Substrings to match to sample names to identify negative controls. Separated by a , | False         | String  | neg,ntc,blank |                                                  |
| --skip_negative_grading     | Skip grading negative controls and just output a PASS for Run QC                    | False         | Boolean | False         |                                                  |

### Other Settings and Parameter Files

Note that the pipeline will create the following files in your current working directory no matter what:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run phac-nml/measeq -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: "./samplesheet.csv"
outdir: "./results/"
reference: "./REFERENCE.fa"
platform: "illumina"
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull phac-nml/measeq
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [phac-nml/measeq releases page](https://github.com/phac-nml/measeq/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow Arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom Configuration

### Resource Requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

## Running in the Background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow Memory Requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
