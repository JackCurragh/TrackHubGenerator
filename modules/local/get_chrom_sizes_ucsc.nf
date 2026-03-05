process GET_CHROM_SIZES_UCSC {
    tag "$ucsc_genome_db"
    label 'process_low'

    container 'docker://quay.io/biocontainers/ucsc-fetchchromsizes:377--h2a80c09_2'

    input:
    val(ucsc_genome_db) // eg "hg38"

    output:
    path("*.chrom_sizes"), emit: chrom_sizes

    script:
    """
    fetchChromSizes ${ucsc_genome_db} > ${ucsc_genome_db}.chrom_sizes
    """

    stub:
    """
    touch ${ucsc_genome_db}.chrom_sizes
    """
}
