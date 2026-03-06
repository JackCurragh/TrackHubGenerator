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
    # Harmonize bedGraph chrom names to match chrom.sizes (ensure UCSC-compat when sizes use chr*)
    awk 'BEGIN{OFS="\t"} FNR==NR {sizes[\$1]=1; next} { if(\$0 ~ /^(track|browser|#)/) next; if(!seen[\$1]++){ if(!( \$1 in sizes)) miss[\$1]=1 }} END{for(c in miss) print c}' \
        $chrom_sizes $bedgraph > ${prefix}.missing.pre || true

    if [ ! -s ${prefix}.missing.pre ]; then
        cp $bedgraph $bgHarmonized
    else
        if awk 'NR==1 { exit (\$1 ~ /^chr/ ? 0 : 1) }' $chrom_sizes; then
            awk 'BEGIN{OFS="\t"} { if(\$0 ~ /^(track|browser|#)/){ print; next } c=\$1; if(c !~ /^chr/){ if(c=="MT"||c=="M"){c="chrM"} else {c="chr" c} } \$1=c; print }' \
                $bedgraph > $bgHarmonized
        else
            cp $bedgraph $bgHarmonized
        fi

        awk 'BEGIN{OFS="\t"} FNR==NR {sizes[\$1]=1; next} { if(\$0 ~ /^(track|browser|#)/) next; if(!seen[\$1]++){ if(!( \$1 in sizes)) miss[\$1]=1 }} END{for(c in miss) print c}' \
            $chrom_sizes $bgHarmonized > ${prefix}.missing.post || true
        if [ -s ${prefix}.missing.post ]; then
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
