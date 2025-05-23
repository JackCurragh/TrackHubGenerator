

process VALIDATE_BEDGRAPH {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::package_name=X.Y.Z"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/package_name:X.Y.Z--hash' :
        'biocontainers/package_name:X.Y.Z--hash' }"

    input:
    tuple val(meta), path(input_file)
    path reference_file
    val(optional_parameter)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Command goes here
    your_command \\
        $args \\
        -i $input_file \\
        -r $reference_file \\
        -o ${prefix}.output

    # Version reporting
    echo '"${task.process}":' > versions.yml
    echo ' command: "$(your_command --version 2>&1 | sed 's/^/    /')"' >> versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bai
    echo '"${task.process}":' > versions.yml
    echo ' command: "your_command_stub"' >> versions.yml
    """
}