process GET_CHROM_SIZES_FASTA {
    tag "$fasta.baseName"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container 'biocontainers/samtools:0.1.19--h94a8ba4_5'

    input:
    path fasta

    output:
    path("chrom.sizes")         , emit: chrom_sizes
    path "versions.yml"         , emit: versions

    script:
    """
    samtools faidx ${fasta} 
    cut -f1,2 ${fasta}.fai > chrom.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch chrom.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
