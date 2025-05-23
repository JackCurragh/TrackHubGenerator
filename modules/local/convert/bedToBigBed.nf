

process CONVERT_BED_TO_BIGBED {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ucsc-bedtobigbed"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedtobigbed:357--1' :
        'biocontainers/ucsc-bedtobigbed:357--1' }"

    input:
    tuple val(meta), path(bed)
    path chrom_sizes

    output:
    tuple val(meta), path("*.bb"), emit: bigwig
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedToBigBed \\
        $args \\
        $bed \\
        $chrom_sizes \\
        ${prefix}.bb

    # Version reporting
    echo '"${task.process}":' > versions.yml
    echo ' command: "$(bedToBigBed --version 2>&1 | sed 's/^/    /')"' >> versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw
    echo '"${task.process}":' > versions.yml
    echo ' command: "your_command_stub"' >> versions.yml
    """
}