process GET_CHROM_SIZES_ASSEMBLY_REPORT {
    tag "$report.baseName"
    label 'process_low'

    container "${ params.container_python ?: 'python:3.11-slim' }"

    input:
    path report
    val  gff_hint

    output:
    path("chrom.sizes") , emit: chrom_sizes
    path "versions.yml"  , emit: versions

    script:
    """
    python - <<'PY'
import os, sys, json

report = "${report}"
gff_hint = r"""${gff_hint}"""

def sniff_chr_from_gff(path):
    if not path or not os.path.exists(path):
        return False
    import gzip
    import io
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            seqid = line.split('\t', 1)[0].strip()
            return seqid.lower().startswith('chr')
    return False

want_ucsc = sniff_chr_from_gff(gff_hint)

headers = []
with open(report, 'rt', encoding='utf-8', errors='ignore') as fh:
    for line in fh:
        if line.startswith('#'):
            # Track the last header-like comment that lists columns
            if 'Sequence-Name' in line and ('Sequence-Length' in line or 'sequence length' in line.lower()):
                # Normalize separators to tabs
                h = line.lstrip('#').strip()
                h = '\t'.join(h.split()) if '\t' not in h else h
                headers = [x.strip() for x in h.split('\t')]

# Build index mapping with normalized keys
idx = {}
for i, h in enumerate(headers):
    key = h.strip().lower().replace(' ', '-').replace('_', '-')
    idx[key] = i

seqname_key = 'sequence-name'
len_key = 'sequence-length'
ucsc_key = 'ucsc-style-name'

use_ucsc = want_ucsc and ucsc_key in idx

def maybe_chrize(name):
    # Add chr prefix for canonical chromosomes if needed
    n = name
    if n in {'X','Y','M','MT'} or n.isdigit():
        return 'chr' + ('M' if n in {'M','MT'} else n)
    return n

out = open('chrom.sizes', 'w')
with open(report, 'rt', encoding='utf-8', errors='ignore') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 2:
            continue
        # Fallback to positional indices if header not parsed
        if idx:
            name = parts[idx.get(seqname_key, 0)]
            length = parts[idx.get(len_key, -1)]
            ucsc = parts[idx.get(ucsc_key, 0)] if ucsc_key in idx else 'na'
        else:
            # Heuristic: sequence name in column 1, length near the end
            name = parts[0]
            length = parts[-2] if parts[-2].isdigit() else parts[-1]
            ucsc = 'na'
        if use_ucsc and ucsc != 'na' and ucsc:
            emit_name = ucsc
        elif want_ucsc:
            emit_name = maybe_chrize(name)
        else:
            emit_name = name
        # Ensure length is integer-like
        try:
            l = int(length)
        except Exception:
            continue
        out.write(f"{emit_name}\t{l}\n")
out.close()

open('versions.yml','w').write('"GET_CHROM_SIZES_ASSEMBLY_REPORT":\n  python: "' + sys.version.split()[0] + '"\n')
PY
    """

    stub:
    """
    echo -e "chr1\t248956422" > chrom.sizes
    printf 'stub' > versions.yml
    """
}

