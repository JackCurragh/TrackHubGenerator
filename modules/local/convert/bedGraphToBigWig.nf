

process CONVERT_BEDGRAPH_TO_BIGWIG {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ucsc-bedgraphtobigwig"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig:357--1' :
        'biocontainers/ucsc-bedgraphtobigwig:357--1' }"

    input:
    tuple val(meta), path(bedgraph)
    path chrom_sizes

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedGraphToBigWig \\
        $args \\
        $bedgraph \\
        $chrom_sizes \\
        ${prefix}.bw

    # Version reporting
    echo '"${task.process}":' > versions.yml
    echo ' command: "$(bedGraphToBigWig --version 2>&1 | sed 's/^/    /')"' >> versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw
    echo '"${task.process}":' > versions.yml
    echo ' command: "your_command_stub"' >> versions.yml
    """
}