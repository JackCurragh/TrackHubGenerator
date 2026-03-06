process CONVERT_BEDGRAPH_TO_BIGWIG {
    tag "$meta.id"
    label 'process_low'

    container "docker://quay.io/biocontainers/ucsc-bedgraphtobigwig:357--1"

    input:
    tuple val(meta), path(bedgraph), path(chrom_sizes)

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bgHarmonized = "${prefix}.harmonized.bedgraph"
    def bgSorted = "${prefix}.sorted.bedgraph"
    """
    # Harmonize bedGraph chrom names to match chrom.sizes per-contig (only change when it yields a match)
    awk '
        BEGIN{OFS="\t"}
        FNR==NR {sizes[\$1]=1; next}
        {
            if(\$0 ~ /^(track|browser|#)/){ print; next }
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
        $chrom_sizes $bedgraph > $bgHarmonized

    # Validate or optionally drop missing contigs (ignore header lines)
    awk 'BEGIN{OFS="\t"} FNR==NR {sizes[\$1]=1; next} { if(\$0 ~ /^(track|browser|#)/) next; if(!seen[\$1]++){ if(!( \$1 in sizes)) miss[\$1]=1 }} END{for(c in miss) print c}' \
        $chrom_sizes $bgHarmonized > ${prefix}.missing.post || true
    if [ -s ${prefix}.missing.post ]; then
        if ${params.drop_missing_contigs ?: true}; then
            echo "[WARN] Dropping records on contigs absent from chrom.sizes:" >&2
            head -n 50 ${prefix}.missing.post >&2
            awk 'BEGIN{OFS="\t"} FNR==NR {sizes[\$1]=1; next} (NR==FNR){next} { if(\$0 ~ /^(track|browser|#)/){ print; next } if(\$1 in sizes) print }' $chrom_sizes $bgHarmonized > ${bgHarmonized}.filtered
            mv ${bgHarmonized}.filtered $bgHarmonized
        else
            echo "[ERROR] The following contigs in bedGraph are absent from chrom.sizes:" >&2
            head -n 50 ${prefix}.missing.post >&2
            exit 2
        fi
    fi

    # Sort bedGraph (required by bedGraphToBigWig); drop headers/comments
    awk 'BEGIN{OFS="\t"} !/^(track|browser|#)/ {print}' $bgHarmonized | sort -k1,1 -k2,2n > $bgSorted

    bedGraphToBigWig \
        $args \
        $bgSorted \
        $chrom_sizes \
        ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedgraphtobigwig: \$(bedGraphToBigWig 2>&1 | head -n1 | sed 's/^bedGraphToBigWig v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedgraphtobigwig: \$(bedGraphToBigWig 2>&1 | head -n1 | sed 's/^bedGraphToBigWig v//')
    END_VERSIONS
    """
}
