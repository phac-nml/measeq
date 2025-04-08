# MeaSeq pipeline: Usage

## Introduction
This pipeline is intended to be run on measles virus paired-end and single-end illumina sequencing data. It customly configures the [viralrecon pipeline](https://nf-co.re/viralrecon/2.6.0) to process measles virus data and follows it with various custom subanalyses that output final variant metrics and highlight mutatations. The pipeline is intended for rapid deployment in outbreak situation in Canada and abroad.

## Index
## Index
- [Running the Pipeline](#running-the-pipeline)
    - [Pipeline Modes](#pipeline-modes)
    - [Mode: `full`](#mode-full)
        - [Profiles](#profiles)
        - [Additional Nextflow Config File](#additional-nextflow-config-file)
        - [Mode: `full` with preset strain reference](#mode-full-with-preset-strain-reference)
            - [Required Input Parameters](#required-input-parameters)
        - [Mode: `full` with user supplied reference and annotation file](#mode-full-with-user-supplied-reference-and-annotation-file)
            - [Required Input Parameters](#required-input-parameters-1)
    - [Mode: `report`](#mode-report)
        - [Required Input Parameters](#required-input-parameters-2)
    - [Optional Parameters](#optional-parameters-for-both-modes)
    - [Fastq File Formatting](#fastq-file-formatting)
    - [Contact Information](#contact-information)
- [Updating the Pipeline](#updating-the-pipeline)

## Running the Pipeline

### Pipeline Modes

There are two modes when running the pipeline
- `full` (default)
    - This mode runs the full pipeline, including viralrecon. This is the default option when running the pipeline and doesn't need to be specified as an argument.
- `report`
    - This mode runs the pipeline using an existing viralrecon results directory. It needs to be specified when running the pipeline.

---

### Mode: `full`

#### Profiles

Profiles are used to specify dependency installation, resources, and how to handle pipeline jobs. They can be passed with `--profile <PROFILE>` or `-p <PROFILE>` and are required when running pipeline in [full mode](#mode-full-with-predefined-strain) but not required when running in [report mode](#mode-report). Please only use one profile for each run.

> **Note**: Dependency management done with the `conda` or `mamba` profiles will take longer to resolve when running viralrecon

Available:

- `conda`: Utilize conda to install dependencies and environment management
- `mamba`: Utilize mamba to install dependencies and environment management
- `singularity`: Utilize singularity for dependencies and environment management
- `docker`: Utilize docker to for dependencies and environment management

> There are many [nextflow environment variables](https://www.nextflow.io/docs/latest/reference/env-vars.html) that you can access to help streamline the analysis on your system. You can set them with:
> ```bash
> export <ENVIRONMENT_VARIABLE>="Val"
> ```
> Of note, when running with conda, it is recommended to set the `NXF_CONDA_CACHEDIR` variable to a shared location to allow easier environment reuse
---

#### Additional Nextflow Config File

While the pipeline uses a [custom viralrecon configuration](../measeq/src/measeq/scripts/measles_parameters.config) file adjusting some tool parameters and versions, you are able to pass further configurations to the viralrecon run by using `--config <CONFIG_FILE>`.

Example:
```bash
measeq -p <PROFILE> --config <CONFIG_FILE> s <B3/D4/D8/H1> -f <PATH/TO/FASTQ>
```

Where in the config you extend with something like:
```
env {
    OPENBLAS_NUM_THREADS = 1
}
executor {
    name = 'slurm'
    queueSize = 50
}
process {
    // Partition to run on
    queue = "MyLargePartition"
}
```

---

#### Mode: `full` with preset strain reference

##### Required Input Parameters
| Parameter | Description | Accepted Inputs |
| - | - | - |
| `-p` / `--profile` | specify the dependency manager to be used | `conda`/`mamba`/`singularity`/`docker` |
| `-s` / `--strain` | Specify which measles strain to use as reference | `B3`/`D4`/`D8`/`H1` |
| `-f` / `--fastq` | Path to directory containing fastq or fastq.gz files | Path to Directory |

- To specify the directory where your fastq files are located, use the `-f` or `--fastq` argument followed by the directory path. Make sure only your fastq files are saved within that directory. Be sure to include what strain you want your reference to be using `-s` or `--strain`.

Example:
```bash
measeq -p <PROFILE> -s <B3/D4/D8/H1> -f <PATH/TO/FASTQ> -o <OUT/DIRECTORY>
```

---

#### Mode: `full` with user supplied reference and annotation file

##### Required Input Parameters
| Parameter | Description | Accepted Inputs |
| - | - | - |
| `-p` / `--profile` | specify the dependency manager to be used | `conda`/`mamba`/`singularity`/`docker` |
| `-r` / `--reference` | Path to FASTA sequence to be used as a reference | Path to File |
| `-g` / -`-gff` | Path to annotation file | Path to File |
| `-f` / `--fastq` | Path to directory containing fastq or fastq.gz files | Path to Directory |

- When running the pipeline with a user supplied reference FASTA and annotation file, the only other argument required is `-f` or `--fastq` followed by the directory path to where the fastq files are located.

Example:
```bash
measeq -p <PROFILE> -r <REF> -g <GFF> -f <PATH/TO/FASTQ>
```

---

### Mode: `report`

##### Required Input Parameters
| Parameter | Description | Accepted Inputs |
| - | - | - |
| `-r` / -`-reference` | Path to FASTA sequence to be used as a reference | Path to File |
| `-vr` / `--viralrecon` | Path to the results directory of viralrecon | Path to Directory |

- The `report` mode uses an existing viralrecon results directory to execute the remaining MeaSeq pipeline steps. It requires the use of the `-vr` or `--viralrecon` parameter followed by the path to the directory. This mode also requires the user to use the `-r` or `--reference` parameter with the same reference FASTA used for the original viralrecon run. 

Example:
```bash
measeq report -r <REF> -vr <PATH/TO/VR/DIR>
```

---

### Optional parameters for both modes
| Parameter | Description | Accepted Inputs |
| - | - | - |
| `-o` / `--output` | Output directory name/path to write results to. Default: './measeq_run_STRAIN_DATE+TIME' | String name or Path to Directory |
| `-c` / `--cpus` | The number of CPUs to be used for downstream parallel processing. Default: 4 | Integer number of CPUs |
| `-h` / `--help` | Displays the help statement and exits. Can be used with the mode parameter for a more specific help statement | None |
| `-v` / `--version` | Displays the version and exits | None |
| `--debug` | Debug mode prints more information to command line | None |

---

### Fastq file formatting
To run this pipeline in `full` mode, make sure your paired reads are saved in the same folder and have the same naming format with only a change of `R1` to `R2` or vice-versa. For example, you can use the following naming format:
```
MEV-001_L001_R1.fastq.gz
MEV-001_L001_R2.fastq.gz
```
If you are using a different naming format, just ensure they match up except for `R1` & `R2`.

> **Important:** Please ensure only fastq or fastq.gz files are saved within your fastq directory. 

---

### Contact Information

By default, the pipeline will display "Contact Information not provided" in the Run Report. If you would like to add any information about the lab/organization running the pipeline, you can use any combination of the following the parameters.

| Parameter | Description | Accepted Inputs | Required? |
| - | - | - | - |
| `--name` | Name of the lab or organizaton running MeaSeq | String name | No |
| `--email` | Email of the lab or organizaton running MeaSeq | Email address | No |
| `--phone` | Phone Number of the lab or organizaton running MeaSeq | Phone Number | No |
| `--website` | Website of the lab or organizaton running MeaSeq | Website | No |

Example:
```bash
measeq \
    -p <PROFILE> \
    -f <PATH/TO/FASTQ> \
    -s <B3/D4/D8/H1> \
    --name <NAME> \
    --email <NAME@ORG.COM> \
    --phone <+123 456 789> \
    --website <WWW.ORG.COM>
```

---

## Updating the pipeline
From time to time, the pipeline's code will be updated. To ensure that you are running the latest version of MeaSeq, it is important you rerun the installation step every time there is a new version of MeaSeq.
