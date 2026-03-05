#!/usr/bin/env python3
import argparse
import gzip
import os
import sys

def sniff_chr_from_gff(path: str) -> bool:
    if not path or not os.path.exists(path):
        return False
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            seqid = line.split('\t', 1)[0].strip()
            return seqid.lower().startswith('chr')
    return False

def parse_args():
    p = argparse.ArgumentParser(description='Build chrom.sizes from NCBI assembly_report.txt with optional chr logic from a GFF hint.')
    p.add_argument('--report', required=True, help='Path to assembly_report.txt')
    p.add_argument('--gff', default='', help='Path to a representative GFF3 (can be .gz) to sniff chr usage')
    p.add_argument('--out', required=True, help='Output chrom.sizes path')
    return p.parse_args()

def main():
    a = parse_args()
    want_ucsc = sniff_chr_from_gff(a.gff)

    headers = []
    with open(a.report, 'rt', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if line.startswith('#'):
                if 'Sequence-Name' in line and ('Sequence-Length' in line or 'sequence length' in line.lower()):
                    h = line.lstrip('#').strip()
                    if '\t' not in h:
                        h = '\t'.join(h.split())
                    headers = [x.strip() for x in h.split('\t')]

    idx = {}
    for i, h in enumerate(headers):
        key = h.strip().lower().replace(' ', '-').replace('_', '-')
        idx[key] = i

    seqname_key = 'sequence-name'
    len_key = 'sequence-length'
    ucsc_key = 'ucsc-style-name'

    use_ucsc = want_ucsc and ucsc_key in idx

    def maybe_chrize(name: str) -> str:
        n = name
        if n in {'X','Y','M','MT'} or n.isdigit():
            return 'chr' + ('M' if n in {'M','MT'} else n)
        return n

    written = 0
    with open(a.report, 'rt', encoding='utf-8', errors='ignore') as fh, open(a.out, 'w') as out:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            if idx:
                name = parts[idx.get(seqname_key, 0)]
                length = parts[idx.get(len_key, -1)]
                ucsc = parts[idx.get(ucsc_key, 0)] if ucsc_key in idx else 'na'
            else:
                name = parts[0]
                length = parts[-2] if len(parts) >= 2 and parts[-2].isdigit() else parts[-1]
                ucsc = 'na'
            if use_ucsc and ucsc != 'na' and ucsc:
                emit_name = ucsc
            elif want_ucsc:
                emit_name = maybe_chrize(name)
            else:
                emit_name = name
            try:
                lval = int(length)
            except Exception:
                continue
            out.write(f"{emit_name}\t{lval}\n")
            written += 1
    if written == 0:
        print('No rows written to chrom.sizes; check assembly_report format', file=sys.stderr)

if __name__ == '__main__':
    main()

