#!/usr/bin/env python3

"""
build_aggregate_manifest.py

Scan a trackhubs output directory and generate a manifest JSON suitable for
the AGGREGATE_TRACKHUB workflow entry. Each manifest entry contains a genome
identifier and the absolute path to its trackDb.txt file.

Example:
  python bin/build_aggregate_manifest.py \
    /path/to/results/trackhubs \
    --output /path/to/results/aggregate_manifest.json

Then run:
  nextflow run main.nf -entry AGGREGATOR \
    --outdir /path/to/results \
    --aggregate_name HPRC_all \
    --aggregate_manifest /path/to/results/aggregate_manifest.json \
    --email you@example.com
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from typing import List, Dict


def find_trackdb_entries(root: str) -> List[Dict[str, str]]:
    """Return a list of {genome, trackdb} entries under the given root.

    This looks for files matching */*/trackDb.txt by default (i.e. one level
    of hub directory then a genome directory), and falls back to a recursive
    search if nothing is found.
    """
    root = os.path.abspath(root)

    # Prefer shallow search (expected layout: <root>/<hub>/<genome>/trackDb.txt)
    matches = glob.glob(os.path.join(root, '*', '*', 'trackDb.txt'))
    if not matches:
        # Fallback: recursive search in case layout differs
        matches = glob.glob(os.path.join(root, '**', 'trackDb.txt'), recursive=True)

    entries: List[Dict[str, str]] = []
    for tdb in matches:
        genome = os.path.basename(os.path.dirname(tdb))
        entries.append({
            'genome': genome,
            'trackdb': os.path.abspath(tdb)
        })

    # De-duplicate and sort for stability
    seen = set()
    deduped: List[Dict[str, str]] = []
    for e in entries:
        key = (e['genome'], e['trackdb'])
        if key not in seen:
            seen.add(key)
            deduped.append(e)
    deduped.sort(key=lambda x: (x['genome'], x['trackdb']))
    return deduped


def main() -> None:
    ap = argparse.ArgumentParser(
        description='Generate aggregate manifest JSON from per-assembly track hubs'
    )
    ap.add_argument('root', help='Root directory containing per-assembly hubs (the "trackhubs" folder)')
    ap.add_argument('-o', '--output', help='Output JSON path (default: <root>/../aggregate_manifest.json)')
    ap.add_argument('-n', '--dry-run', action='store_true', help='Print summary but do not write a file')
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    if not os.path.isdir(root):
        ap.error(f'Root directory does not exist: {root}')

    entries = find_trackdb_entries(root)
    if not entries:
        ap.error(f'No trackDb.txt files found under {root}')

    print(f'Found {len(entries)} genomes under {root}')

    if args.dry_run:
        # Show a small sample
        for e in entries[:5]:
            print(f"  - {e['genome']}: {e['trackdb']}")
        return

    out_default = os.path.abspath(os.path.join(root, os.pardir, 'aggregate_manifest.json'))
    out_path = os.path.abspath(args.output) if args.output else out_default
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with open(out_path, 'w') as fh:
        json.dump(entries, fh, indent=2)

    print(f'Wrote {out_path} with {len(entries)} entries')
    print('Next: use --aggregate_manifest with the AGGREGATOR entry in main.nf')


if __name__ == '__main__':
    main()

