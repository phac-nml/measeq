// Process to just stage the predictions file so that the IRIDA Next plugin can get it
process STAGE_FILE_IRIDANEXT {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/coreutils:8.31--h14c3975_0'
        : 'biocontainers/coreutils:8.31--h14c3975_0' }"

    input:
    path csv

    output:
    path "predictions.csv", includeInputs: true

    script:
    """
    """

    stub:
    """
    """
}
