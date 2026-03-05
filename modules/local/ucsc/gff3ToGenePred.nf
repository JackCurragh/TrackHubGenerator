process UCSC_GFF3_TO_GENEPRED {
    tag "$meta.id"
    label 'gff3_heavy'

    container "docker://quay.io/biocontainers/ucsc-gff3togenepred:377--h2a80c09_2"

    input:
    tuple val(meta), path(gff3)

    output:
    tuple val(meta), path("*.genePred"), emit: genepred

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    IN="${gff3}"
    if [[ "${gff3}" == *.gz ]]; then
        echo "Decompressing ${gff3} -> ${prefix}.gff3" >&2
        gunzip -c "${gff3}" > "${prefix}.gff3"
        IN="${prefix}.gff3"
    fi

    # Strip noisy attribute to reduce warnings, then convert
    cat "\$IN" | awk -F'\t' 'BEGIN{OFS="\t"} /^#/ {print; next} { sub(/;AdditionalSources=[^;]*/, "", $9); print }' > "${prefix}.sanitized.gff3"
    gff3ToGenePred -useName -warnAndContinue "${prefix}.sanitized.gff3" "${prefix}.genePred"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genePred
    """
}
