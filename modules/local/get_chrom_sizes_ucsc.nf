process GET_CHROM_SIZES_UCSC {
    tag "$ucsc_genome_db"
    label 'process_low'

    conda "bioconda::ucsc-fetchchromsizes"
    container 'biocontainers/ucsc-fetchchromsizes:357--1'

    input:
    val(ucsc_genome_db) // eg "hg38"

    output:
    path("*.chrom_sizes"), emit: chrom_sizes
    path "versions.yml"  , emit: versions

    script:
    """
    fetchChromSizes ${ucsc_genome_db} > ${ucsc_genome_db}.chrom_sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc-fetchchromsizes: \$(fetchChromSizes 2>&1 | head -n1 | sed 's/^fetchChromSizes v//' || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${ucsc_genome_db}.chrom_sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc-fetchchromsizes: \$(fetchChromSizes 2>&1 | head -n1 | sed 's/^fetchChromSizes v//' || echo "unknown")
    END_VERSIONS
    """
}
