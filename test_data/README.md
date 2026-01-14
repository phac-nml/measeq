# Test Data

Data that can be used to test that the pipeline is functioning on your system/workflow along with for CI tests.

## Illumina Fastq Data

The [illumina_fastqs](./illumina_fastqs/) directory contains two downsampled and properly paired MeV genotype D8 samples from [SRA BioProject PRJNA480551](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA480551) that can be used for testing the MeaSeq pipeline and for CI tests

To run the test data easily, the pipeline can be run with the `test_illumina` profile like:

```bash
nextflow run phac-nml/measeq -profile test_illumina,<docker/singularity>
```

Results end up in the `results` output directory

## Nanopore Fastq Data

The [nanopore_fastqs](./nanopore_fastqs/) directory contains one downsampled MeV genotype D8 samples from [SRA BioProject PRJNA1174053](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1174053) and one downsampled MeV genotype B3 from [SRA BioProject PRJNA843031](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA843031) that can be used for testing the MeaSeq pipeline and for CI tests

To run the test data easily, the pipeline can be run with the `test_nanopore` profile like:

```bash
nextflow run phac-nml/measeq -profile test_nanopore,<docker/singularity>
```

Results end up in the `results` output directory

## Primer Bed File for Amplicon Testing

There is also an example [primer bed](./primer.bed) file to test out running the amplicon steps for both pipelines. The included file is a test file and should not be used for real analyses.

For internal sequencing we are using a custom amplicon scheme (to be released).
