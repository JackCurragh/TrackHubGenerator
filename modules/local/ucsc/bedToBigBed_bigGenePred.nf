process UCSC_BED_TO_BIGBED_BIGGENEPRED {
    tag "$meta.id"
    label 'process_low'

    container "docker://quay.io/biocontainers/ucsc-bedtobigbed:357--1"

    input:
    tuple val(meta), path(bed), path(chrom_sizes), path(as_file)

    output:
    tuple val(meta), path("*.bb"), emit: bigbed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bedSorted = "${prefix}.sorted.bed"
    def bedHarmonized = "${prefix}.harmonized.bed"
    """
    # Harmonize BED chrom names to match chrom.sizes (ensure UCSC-compat when sizes use chr*)
    # 1) If all contigs already present in sizes -> pass through
    # 2) Else, if sizes use chr-prefix and BED does not -> add chr (and map MT/M->chrM)
    # 3) Recheck; if still mismatched, print offending contigs and exit

    awk 'BEGIN{OFS="\t"} FNR==NR {sizes[\$1]=1; next} {if(!seen[\$1]++){ if(!( \$1 in sizes)) miss[\$1]=1}} END{for(c in miss) print c}' \
        ${chrom_sizes} ${bed} > ${prefix}.missing.pre || true

    if [ ! -s ${prefix}.missing.pre ]; then
        cp ${bed} ${bedHarmonized}
    else
        # Do we have chr-prefixed sizes?
        if awk 'NR==1 { exit (\$1 ~ /^chr/ ? 0 : 1) }' ${chrom_sizes}; then
            awk 'BEGIN{OFS="\t"} {c=$1; if(c !~ /^chr/){ if(c=="MT"||c=="M"){c="chrM"} else {c="chr" c} } $1=c; print}' \
                ${bed} > ${bedHarmonized}
        else
            # Sizes are not chr-prefixed; leave as-is (we do not auto-strip chr here)
            cp ${bed} ${bedHarmonized}
        fi

        # Recheck after attempted harmonization
        awk 'BEGIN{OFS="\t"} FNR==NR {sizes[\$1]=1; next} {if(!seen[\$1]++){ if(!( \$1 in sizes)) miss[\$1]=1}} END{for(c in miss) print c}' \
            ${chrom_sizes} ${bedHarmonized} > ${prefix}.missing.post || true

        if [ -s ${prefix}.missing.post ]; then
            echo "[ERROR] The following contigs in BED are absent from chrom.sizes:" >&2
            head -n 50 ${prefix}.missing.post >&2
            echo "Hint: Provide an AssemblyReport in the samplesheet so chrom.sizes can be built with UCSC-style names, or ensure GFF/BED chrom names match." >&2
            exit 2
        fi
    fi

    sort -k1,1 -k2,2n ${bedHarmonized} > ${bedSorted}
    bedToBigBed -type=bed12+8 -as=${as_file} ${bedSorted} ${chrom_sizes} ${prefix}.bb

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bb

    """
}
