process UCSC_GENEPRED_TO_BIGGENEPRED_BED {
    tag "$meta.id"
    label 'process_low'

    container "biocontainers/ucsc-genepredtobiggenepred:377--h199ee4e_0"

    input:
    tuple val(meta), path(genepred)

    output:
    tuple val(meta), path("*.bigGenePred.bed"), emit: biggenepred_bed
    path "versions.yml"                        , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    genePredToBigGenePred ${genepred} ${prefix}.bigGenePred.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genepredtobiggenepred: \$(genePredToBigGenePred 2>&1 | head -n1 | sed 's/^genePredToBigGenePred v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bigGenePred.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genepredtobiggenepred: "stub"
    END_VERSIONS
    """
}
