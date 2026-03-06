process CONVERT_BED_TO_BIGBED {
    tag "$meta.id"
    label 'process_low'

    container "docker://quay.io/biocontainers/ucsc-bedtobigbed:357--1"

    input:
    tuple val(meta), path(bed), path(chrom_sizes)

    output:
    tuple val(meta), path("*.bb"), emit: bigbed
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bedHarmonized = "${prefix}.harmonized.bed"
    def bedSorted = "${prefix}.sorted.bed"
    """
    # Harmonize BED chrom names to match chrom.sizes per-contig (only change when it yields a match)
    awk '
        BEGIN{OFS="\t"}
        FNR==NR {sizes[\$1]=1; next}
        {
            c=\$1; mapped=c
            if(!(c in sizes)) {
                add = (c=="MT"||c=="M")?"chrM":"chr" c
                if(add in sizes) mapped=add
                else if(c ~ /^chr/) {
                    base=substr(c,4); if(base=="M"||base=="MT") base="MT"
                    if(base in sizes) mapped=base
                }
            }
            \$1=mapped; print
        }' \
        $chrom_sizes $bed > $bedHarmonized

    # Validate or optionally drop missing contigs
    awk 'BEGIN{OFS="\t"} FNR==NR {sizes[\$1]=1; next} {if(!seen[\$1]++){ if(!( \$1 in sizes)) miss[\$1]=1}} END{for(c in miss) print c}' \
        $chrom_sizes $bedHarmonized > ${prefix}.missing.post || true
    if [ -s ${prefix}.missing.post ]; then
        if ${params.drop_missing_contigs ?: true}; then
            echo "[WARN] Dropping records on contigs absent from chrom.sizes:" >&2
            head -n 50 ${prefix}.missing.post >&2
            awk 'BEGIN{OFS="\t"} FNR==NR {sizes[\$1]=1; next} (\$1 in sizes)' $chrom_sizes $bedHarmonized > ${bedHarmonized}.filtered
            mv ${bedHarmonized}.filtered $bedHarmonized
        else
            echo "[ERROR] The following contigs in BED are absent from chrom.sizes:" >&2
            head -n 50 ${prefix}.missing.post >&2
            exit 2
        fi
    fi

    sort -k1,1 -k2,2n $bedHarmonized > $bedSorted

    bedToBigBed \
        $args \
        $bedSorted \
        $chrom_sizes \
        ${prefix}.bb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtobigbed: \$(bedToBigBed 2>&1 | head -n1 | sed 's/^bedToBigBed v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtobigbed: \$(bedToBigBed 2>&1 | head -n1 | sed 's/^bedToBigBed v//')
    END_VERSIONS
    """
}
