process GET_CHROM_SIZES_UCSC {
    tag "$ucsc_genome_db"
    label 'process_low'

    conda "bioconda::ucsc-fetchchromsizes"
    container 'biocontainers/ucsc-fetchchromsizes:357--1'

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
