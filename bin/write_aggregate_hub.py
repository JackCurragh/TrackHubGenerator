#!/usr/bin/env python3

"""
write_aggregate_hub.py

Create a parent (aggregate) track hub from a manifest JSON that lists
per-assembly trackDb.txt files. Intended to be called from Nextflow.

Manifest format: a JSON array of objects with keys:
  - genome: genome identifier (e.g., GCA_...)
  - trackdb: absolute path to that assembly's trackDb.txt

This script writes:
  ./<aggregate_name>/hub.txt
  ./<aggregate_name>/genomes.txt

The genomes.txt file uses relative paths to each trackDb.txt resolving from
<outdir>/trackhubs/<aggregate_name>/ to ensure the hub remains relocatable
under the same root.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from typing import Any, Dict, List


def _norm_opt(value: str | None) -> str | None:
    if value is None:
        return None
    v = str(value).strip()
    if not v or v.lower() == 'null':
        return None
    return v


def _sample_label_from_trackdb(trackdb_abs: str, genome: str) -> str:
    hub_dir = os.path.basename(os.path.dirname(os.path.dirname(trackdb_abs)))
    match = re.match(r'^(?:HPRC_?)?(?P<sample>[A-Z0-9]+)_(?P<hap>mat|pat|hap1|hap2)$', hub_dir)
    if not match:
        return genome
    hap = match.group('hap')
    hap_label = {'mat': 'maternal', 'pat': 'paternal', 'hap1': 'hap1', 'hap2': 'hap2'}[hap]
    return f"{match.group('sample')} {hap_label} ({genome})"


def main() -> None:
    ap = argparse.ArgumentParser(description='Write aggregate parent hub from a manifest JSON')
    ap.add_argument('--manifest', required=True, help='Path to manifest JSON with {genome,trackdb} entries')
    ap.add_argument('--aggregate-name', required=True, help='Aggregate hub name')
    ap.add_argument('--outdir', required=True, help='Pipeline outdir (parent of trackhubs)')
    ap.add_argument('--email', required=True, help='Contact email for hub.txt')
    ap.add_argument('--short-label', default='', help='Short label (<=17 chars); defaults to aggregate name')
    ap.add_argument('--long-label', default='', help='Long label; defaults to "<aggregate-name> multi-assembly hub"')
    args = ap.parse_args()

    agg_name = args.aggregate_name
    outdir = os.path.abspath(args.outdir)
    manifest_path = os.path.abspath(args.manifest)

    with open(manifest_path, 'r') as fh:
        hubs: List[Dict[str, Any]] = json.load(fh)

    short = _norm_opt(args.short_label) or agg_name[:17]
    long = _norm_opt(args.long_label) or f"{agg_name} multi-assembly hub"

    # Where the parent hub directory will live once published
    root = os.path.join(outdir, 'trackhubs')
    target = os.path.join(root, agg_name)

    # Create output directory in the current working directory (Nextflow publishes it)
    agg_dir = os.path.join(agg_name)
    os.makedirs(agg_dir, exist_ok=True)

    # hub.txt
    with open(os.path.join(agg_dir, 'hub.txt'), 'w') as fh:
        fh.write(f"hub {agg_name}\n")
        fh.write(f"shortLabel {short}\n")
        fh.write(f"longLabel {long}\n")
        fh.write("genomesFile genomes.txt\n")
        fh.write(f"email {args.email}\n")
        fh.write("descriptionUrl hubDescription.html\n")

    # genomes.txt
    with open(os.path.join(agg_dir, 'genomes.txt'), 'w') as fh:
        for order, entry in enumerate(hubs, start=1):
            genome = entry['genome']
            trackdb_abs = os.path.abspath(entry['trackdb'])
            rel = os.path.relpath(trackdb_abs, start=target)
            description = _norm_opt(entry.get('description')) or _sample_label_from_trackdb(trackdb_abs, genome)
            fh.write(f"genome {genome}\n")
            fh.write("organism Human\n")
            fh.write("scientificName Homo sapiens\n")
            fh.write(f"description {description}\n")
            fh.write(f"orderKey {order}\n")
            fh.write(f"trackDb {rel}\n\n")

    with open(os.path.join(agg_dir, 'hubDescription.html'), 'w') as fh:
        fh.write("<!doctype html>\n")
        fh.write("<html><head><meta charset=\"utf-8\"><title>HPRC track hub</title></head><body>\n")
        fh.write(f"<h1>{long}</h1>\n")
        fh.write("<p>HPRC annotation tracks for GenArk assemblies. Assemblies are identified by their GCA accessions and labelled by sample and haplotype where this can be inferred from the hub path.</p>\n")
        fh.write("<p>Tracks are generated from CAT and Ensembl GFF3 annotations converted to BigBed/bigGenePred.</p>\n")
        fh.write("</body></html>\n")

    # Optional: small versions file for traceability
    try:
        import platform, sys
        with open('versions.yml', 'w') as vf:
            vf.write('write_aggregate_hub:\n')
            vf.write(f"  python: {sys.version.split()[0]}\n")
            vf.write(f"  platform: {platform.platform()}\n")
    except Exception:
        pass


if __name__ == '__main__':
    main()
