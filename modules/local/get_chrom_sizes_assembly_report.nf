process GET_CHROM_SIZES_ASSEMBLY_REPORT {
    tag "$genome"
    label 'process_low'

    container "python:3.11-slim"

    input:
    tuple val(genome), path(report), val(gff_hint)

    output:
    tuple val(genome), path("chrom.sizes"), emit: chrom_sizes

    script:
    def gff_arg = (gff_hint && gff_hint.toString().trim()) ? "--gff '${gff_hint}'" : ''
    def force_ucsc = params.force_ucsc_chrom_names ? '--force-ucsc' : ''
    """
    assembly_report_to_chrom_sizes.py \
      --report '${report}' \
      ${gff_arg} \
      ${force_ucsc} \
      --out chrom.sizes
    """

    stub:
    """
    echo -e "chr1\t248956422" > chrom.sizes
    printf 'stub' > versions.yml
    """
}
