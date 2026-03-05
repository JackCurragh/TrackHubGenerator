process GET_CHROM_SIZES_FASTA {
    tag "$fasta.baseName"
    label 'process_low'

    container 'docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    input:
    path fasta

    output:
    path("chrom.sizes")         , emit: chrom_sizes

    script:
    """
    samtools faidx ${fasta} 
    cut -f1,2 ${fasta}.fai > chrom.sizes
    """

    stub:
    """
    touch chrom.sizes
    """
}
