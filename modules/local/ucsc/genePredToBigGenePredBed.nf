process UCSC_GENEPRED_TO_BIGGENEPRED_BED {
    tag "$meta.id"
    label 'process_low'

    container "docker://quay.io/biocontainers/ucsc-genepredtobiggenepred:377--h2a80c09_2"

    input:
    tuple val(meta), path(genepred)

    output:
    tuple val(meta), path("*.bigGenePred.bed"), emit: biggenepred_bed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    genePredToBigGenePred ${genepred} ${prefix}.bigGenePred.bed

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bigGenePred.bed

    """
}
