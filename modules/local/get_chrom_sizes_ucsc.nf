

process PROCESS_NAME {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ucsc-fetchchromsizes"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-fetchchromsizes:357--1' :
        'biocontainers/ucsc-fetchchromsizes:357--1' }"

    input:
    val(ucsc_genome_db) // eg "hg38"

    output:
    path("*.chrom_sizes"), emit: bam
    path "versions.yml"  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    fetchChromSizes ${ucsc_genome_db} > ${ucsc_genome_db}.chrom_sizes

    # Version reporting
    echo '"${task.process}":' > versions.yml
    echo ' command: "$(fetchChromSizes --version 2>&1 | sed 's/^/    /')"' >> versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${ucsc_genome_db}.chrom_sizes
    echo '"${task.process}":' > versions.yml
    echo ' command: "your_command_stub"' >> versions.yml
    """
}