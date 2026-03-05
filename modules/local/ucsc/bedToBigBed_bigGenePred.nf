process UCSC_BED_TO_BIGBED_BIGGENEPRED {
    tag "$meta.id"
    label 'process_low'

    container "biocontainers/ucsc-bedtobigbed:357--1"

    input:
    tuple val(meta), path(bed), path(chrom_sizes), path(as_file)

    output:
    tuple val(meta), path("*.bb"), emit: bigbed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bedSorted = "${prefix}.sorted.bed"
    """
    sort -k1,1 -k2,2n ${bed} > ${bedSorted}
    bedToBigBed -type=bed12+8 -as=${as_file} ${bedSorted} ${chrom_sizes} ${prefix}.bb

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bb

    """
}
