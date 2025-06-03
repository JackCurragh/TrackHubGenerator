process CONVERT_BEDGRAPH_TO_BIGWIG {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ucsc-bedgraphtobigwig"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/ucsc-bedgraphtobigwig:481--6675980cef0e7276' :
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedgraphtobigwig: \$(bedGraphToBigWig 2>&1 | head -n1 | sed 's/^bedGraphToBigWig v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedgraphtobigwig: \$(bedGraphToBigWig 2>&1 | head -n1 | sed 's/^bedGraphToBigWig v//')
    END_VERSIONS
    """
}