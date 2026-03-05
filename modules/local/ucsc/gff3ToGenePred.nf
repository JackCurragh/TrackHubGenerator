process UCSC_GFF3_TO_GENEPRED {
    tag "$meta.id"
    label 'process_low'

    container "${ params.container_ucsc_gff3togenepred ?: 'biocontainers/ucsc-gff3togenepred' }"

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("*.genePred"), emit: genepred
    path "versions.yml"                 , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    IN=${gff3}
    if [[ "${gff3}" == *.gz ]]; then
        echo "Decompressing ${gff3} -> ${prefix}.gff3" >&2
        gunzip -c ${gff3} > ${prefix}.gff3
        IN=${prefix}.gff3
    fi

    gff3ToGenePred -useName -warnAndContinue ${IN} ${prefix}.genePred

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gff3togenepred: $(gff3ToGenePred 2>&1 | head -n1 | sed 's/^gff3ToGenePred v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genePred
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gff3togenepred: "stub"
    END_VERSIONS
    """
}
