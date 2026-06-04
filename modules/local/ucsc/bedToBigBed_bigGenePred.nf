process UCSC_BED_TO_BIGBED_BIGGENEPRED {
    tag "$meta.id"
    label 'process_low'

    container "quay.io/biocontainers/ucsc-bedtobigbed:357--1"

    input:
    tuple val(meta), path(bed), path(chrom_sizes), path(as_file), path(assembly_report)

    output:
    tuple val(meta), path("*.bb"), emit: bigbed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bedSorted = "${prefix}.sorted.bed"
    def bedHarmonized = "${prefix}.harmonized.bed"
    def bedRemapped = "${prefix}.remapped.bed"
    def chromColumn = params.chrom_name_column ?: 'auto'
    def dropMissing = params.drop_missing_contigs == true
    def targetColumnIndex = [
        'sequence-name': 1,
        'assigned-molecule': 3,
        'genbank-accn': 5,
        'refseq-accn': 7,
        'ucsc-style-name': 10
    ][chromColumn] ?: 0
    """
    cp ${bed} ${bedRemapped}

    if [ -s ${assembly_report} ] && [ "${chromColumn}" != "auto" ] && [ "${targetColumnIndex}" != "0" ]; then
        awk -v target_col="${targetColumnIndex}" '
            BEGIN{FS=OFS="\t"}
            FNR==NR {
                if(\$0 ~ /^#/) next
                target=\$target_col
                if(target == "" || target == "na") next
                alias_cols[1]=1; alias_cols[2]=3; alias_cols[3]=5; alias_cols[4]=7; alias_cols[5]=10
                for(i=1; i<=5; i++) {
                    alias=\$alias_cols[i]
                    if(alias != "" && alias != "na") map[alias]=target
                }
                if(\$3 == "MT") map["M"]=target
                if(\$3 == "M") map["MT"]=target
                next
            }
            {
                if(\$1 in map) \$1=map[\$1]
                print
            }
        ' ${assembly_report} ${bed} > ${bedRemapped}
    fi

    # Harmonize BED chrom names to match chrom.sizes per-contig (only change when it yields a match)
    awk '
        BEGIN{FS=OFS="\t"}
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
        ${chrom_sizes} ${bedRemapped} > ${bedHarmonized}

    # Validate or optionally drop missing contigs
    awk 'BEGIN{FS=OFS="\t"} FNR==NR {sizes[\$1]=1; next} {if(!seen[\$1]++){ if(!( \$1 in sizes)) miss[\$1]=1}} END{for(c in miss) print c}' \
        ${chrom_sizes} ${bedHarmonized} > ${prefix}.missing.post || true

    if [ -s ${prefix}.missing.post ]; then
        if ${dropMissing}; then
            echo "[WARN] Dropping records on contigs absent from chrom.sizes:" >&2
            head -n 50 ${prefix}.missing.post >&2
            awk 'BEGIN{FS=OFS="\t"} FNR==NR {sizes[\$1]=1; next} (\$1 in sizes)' ${chrom_sizes} ${bedHarmonized} > ${bedHarmonized}.filtered
            mv ${bedHarmonized}.filtered ${bedHarmonized}
        else
            echo "[ERROR] The following contigs in BED are absent from chrom.sizes:" >&2
            head -n 50 ${prefix}.missing.post >&2
            echo "Hint: Ensure your AssemblyReport/FASTA matches the GFF; otherwise those contigs cannot be converted." >&2
            exit 2
        fi
    fi

    if [ ! -s ${bedHarmonized} ]; then
        echo "[ERROR] No BED rows remain for ${prefix}; refusing to create an empty BigBed." >&2
        exit 3
    fi

    sort -k1,1 -k2,2n ${bedHarmonized} > ${bedSorted}
    bedToBigBed -tab -type=bed12+8 -as=${as_file} ${bedSorted} ${chrom_sizes} ${prefix}.bb

    if command -v bigBedInfo >/dev/null 2>&1; then
        item_count=\$(bigBedInfo ${prefix}.bb 2>/dev/null | awk '/^itemCount:/ {print \$2}')
        if [ "\${item_count:-0}" = "0" ]; then
            echo "[ERROR] BigBed ${prefix}.bb has zero items; refusing to publish it." >&2
            exit 4
        fi
    else
        echo "[WARN] bigBedInfo not available; skipped post-conversion itemCount validation." >&2
    fi

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bb

    """
}
