

process PROCESS_NAME {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::package_name=X.Y.Z"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/package_name:X.Y.Z--hash' :
        'biocontainers/package_name:X.Y.Z--hash' }"

    input:
    path fasta

    output:
    path("chrom.sizes")   , emit: chrom_sizes
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools faidx ${fasta} 
    cut -f1,2 i${fasta}.fai > chrom.sizes

    # Version reporting
    echo '"${task.process}":' > versions.yml
    echo ' command: "$(samtools --version 2>&1 | sed 's/^/    /')"' >> versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch chrom.sizes
    echo '"${task.process}":' > versions.yml
    echo ' command: "your_command_stub"' >> versions.yml
    """
}