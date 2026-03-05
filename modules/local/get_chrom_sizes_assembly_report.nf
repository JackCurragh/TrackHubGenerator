process GET_CHROM_SIZES_ASSEMBLY_REPORT {
    tag "$report.baseName"
    label 'process_low'

    container "python:3.11-slim"

    input:
    path report
    val  gff_hint

    output:
    path("chrom.sizes") , emit: chrom_sizes
    path "versions.yml"  , emit: versions

    script:
    def gff_arg = (gff_hint && gff_hint.toString().trim()) ? "--gff '${gff_hint}'" : ''
    """
    assembly_report_to_chrom_sizes.py \
      --report '${report}' \
      ${gff_arg} \
      --out chrom.sizes
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      python: \$(python -V 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "chr1\t248956422" > chrom.sizes
    printf 'stub' > versions.yml
    """
}
