#!/usr/bin/env python3
import argparse
import gzip
import os
import sys
from typing import Dict, List, Tuple, Set

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

def sample_gff_seqids(path: str, limit: int = 10000) -> Set[str]:
    """Collect up to `limit` seqids from a GFF to guide name-column choice."""
    ids: Set[str] = set()
    if not path or not os.path.exists(path):
        return ids
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            seqid = line.split('\t', 1)[0].strip()
            if seqid:
                ids.add(seqid)
                if len(ids) >= limit:
                    break
    return ids

def parse_args():
    p = argparse.ArgumentParser(description='Build chrom.sizes from NCBI assembly_report.txt with optional chr logic from a GFF hint.')
    p.add_argument('--report', required=True, help='Path to assembly_report.txt')
    p.add_argument('--gff', default='', help='Path to a representative GFF3 (can be .gz) to sniff chr usage')
    p.add_argument('--force-ucsc', action='store_true', help='Force UCSC-style names (chr*); ignores GFF sniffing')
    p.add_argument('--out', required=True, help='Output chrom.sizes path')
    return p.parse_args()

def main():
    a = parse_args()
    want_ucsc = a.force_ucsc or sniff_chr_from_gff(a.gff)
    gff_ids = sample_gff_seqids(a.gff, limit=20000)

    # Parse header
    headers: List[str] = []
    with open(a.report, 'rt', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if line.startswith('#'):
                if 'Sequence-Name' in line and ('Sequence-Length' in line or 'sequence length' in line.lower()):
                    h = line.lstrip('#').strip()
                    if '\t' not in h:
                        h = '\t'.join(h.split())
                    headers = [x.strip() for x in h.split('\t')]

    idx: Dict[str, int] = {}
    for i, h in enumerate(headers):
        key = h.strip().lower().replace(' ', '-').replace('_', '-')
        idx[key] = i

    # Candidate name columns, in preference order
    candidates = [
        'ucsc-style-name',
        'sequence-name',
        'refseq-accn',
        'genbank-accn',
    ]
    len_key = 'sequence-length'

    # Build per-column name sets for match scoring
    name_sets: Dict[str, Set[str]] = {k: set() for k in candidates if k in idx}
    if idx:
        with open(a.report, 'rt', encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                for k in list(name_sets.keys()):
                    v = parts[idx[k]] if idx[k] < len(parts) else ''
                    if v and v.lower() != 'na':
                        name_sets[k].add(v)

    # Choose the best-matching column
    chosen = None
    if name_sets:
        # Score by intersection with GFF ids
        scores = {k: (len(name_sets[k] & gff_ids) if gff_ids else 0) for k in name_sets}
        # If forcing UCSC and ucsc column exists and has any non-na values, bias strongly
        if a.force_ucsc and 'ucsc-style-name' in name_sets and len(name_sets['ucsc-style-name']) > 0:
            # Use UCSC if it matches OR if no other column matches at all
            if (gff_ids and scores['ucsc-style-name'] >= max(scores.values())) or (not gff_ids and len(name_sets['ucsc-style-name']) > 0):
                chosen = 'ucsc-style-name'
        if not chosen:
            # Otherwise pick the highest score; if all zero, fall back to sequence-name, then refseq, then genbank
            if scores:
                best = max(scores, key=lambda k: (scores[k], 1 if k == 'sequence-name' else 0))
                chosen = best
    # Final fallback when no headers/candidates were parsed
    if not chosen:
        chosen = 'sequence-name' if 'sequence-name' in idx else ( 'refseq-accn' if 'refseq-accn' in idx else ( 'genbank-accn' if 'genbank-accn' in idx else None ) )

    def maybe_chrize(name: str) -> str:
        n = name
        if n in {'X','Y','M','MT'} or n.isdigit():
            return 'chr' + ('M' if n in {'M','MT'} else n)
        return n

    written = 0
    seen: Set[str] = set()
    with open(a.report, 'rt', encoding='utf-8', errors='ignore') as fh, open(a.out, 'w') as out:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            if idx and chosen and chosen in idx and len_key in idx and idx[chosen] < len(parts) and idx[len_key] < len(parts):
                raw_name = parts[idx[chosen]]
                if not raw_name or raw_name.lower() == 'na':
                    # fall back to sequence-name if available
                    raw_name = parts[idx.get('sequence-name', idx[chosen])] if idx.get('sequence-name', None) is not None else raw_name
                length = parts[idx[len_key]]
            else:
                raw_name = parts[0]
                length = parts[-2] if len(parts) >= 2 and parts[-2].isdigit() else parts[-1]

            # Apply chr-ization only for simple primary names when UCSC is desired and no explicit UCSC-style exists
            emit_name = raw_name
            if a.force_ucsc and (not chosen or chosen != 'ucsc-style-name'):
                emit_name = maybe_chrize(raw_name)

            try:
                lval = int(length)
            except Exception:
                continue
            if not emit_name or emit_name.lower() == 'na' or emit_name in seen:
                continue
            out.write(f"{emit_name}\t{lval}\n")
            seen.add(emit_name)
            written += 1
    if written == 0:
        print('No rows written to chrom.sizes; check assembly_report format', file=sys.stderr)

if __name__ == '__main__':
    main()
