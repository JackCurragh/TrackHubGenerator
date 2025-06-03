process CONVERT_BED_TO_BIGBED {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::ucsc-bedtobigbed"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/ucsc-bedtobigbed:473--d4d499b685c95583' :
        'biocontainers/ucsc-bedtobigbed:357--1' }"

    input:
    tuple val(meta), path(bed)
    path chrom_sizes

    output:
    tuple val(meta), path("*.bb"), emit: bigbed
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtobigbed: \$(bedToBigBed 2>&1 | head -n1 | sed 's/^bedToBigBed v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtobigbed: \$(bedToBigBed 2>&1 | head -n1 | sed 's/^bedToBigBed v//')
    END_VERSIONS
    """
}