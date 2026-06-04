#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a JBrowse 2 config.json from generated per-assembly track hubs."
    )
    parser.add_argument("--trackhubs-root", required=True, help="Directory containing generated hubs")
    parser.add_argument("--base-url", required=True, help="Public URL corresponding to --trackhubs-root")
    parser.add_argument(
        "--twobit-url-template",
        required=True,
        help="URL template for assembly 2bit files. May contain {genome}.",
    )
    parser.add_argument(
        "--chrom-sizes-url-template",
        default="",
        help="Optional URL template for chrom.sizes files. May contain {genome}.",
    )
    parser.add_argument("--out", required=True, help="Output JBrowse config.json")
    return parser.parse_args()


def parse_trackdb(path: Path) -> list[dict[str, str]]:
    tracks: list[dict[str, str]] = []
    current: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not raw.strip():
            continue
        line = raw.strip()
        parts = line.split(None, 1)
        if len(parts) == 2 and parts[0] == "track":
            if current.get("track") and current.get("bigDataUrl"):
                tracks.append(current)
            current = {}
        if len(parts) == 2:
            current[parts[0]] = parts[1]
    if current.get("track") and current.get("bigDataUrl"):
        tracks.append(current)
    return tracks


def clean_track_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def make_uri(base_url: str, root: Path, file_path: Path) -> str:
    rel = file_path.relative_to(root).as_posix()
    return f"{base_url.rstrip('/')}/{rel}"


def format_template(template: str, genome: str) -> str:
    return template.replace("{genome}", genome)


def main() -> None:
    args = parse_args()
    root = Path(args.trackhubs_root).resolve()
    base_url = args.base_url.rstrip("/")

    assemblies: dict[str, dict[str, Any]] = {}
    tracks: list[dict[str, Any]] = []

    for trackdb in sorted(root.glob("*/*/trackDb.txt")):
        genome = trackdb.parent.name
        if genome not in assemblies:
            adapter: dict[str, Any] = {
                "type": "TwoBitAdapter",
                "twoBitLocation": {
                    "uri": format_template(args.twobit_url_template, genome),
                    "locationType": "UriLocation",
                },
            }
            if args.chrom_sizes_url_template:
                adapter["chromSizesLocation"] = {
                    "uri": format_template(args.chrom_sizes_url_template, genome),
                    "locationType": "UriLocation",
                }
            assemblies[genome] = {
                "name": genome,
                "sequence": {
                    "type": "ReferenceSequenceTrack",
                    "trackId": f"{clean_track_id(genome)}-ReferenceSequenceTrack",
                    "adapter": adapter,
                },
            }

        for stanza in parse_trackdb(trackdb):
            data_url = stanza["bigDataUrl"]
            bigbed_path = (trackdb.parent / data_url).resolve()
            if bigbed_path.exists():
                uri = make_uri(base_url, root, bigbed_path)
            else:
                uri = data_url
            track_id = clean_track_id(f"{genome}-{stanza['track']}")
            tracks.append(
                {
                    "type": "FeatureTrack",
                    "trackId": track_id,
                    "name": stanza.get("longLabel") or stanza.get("shortLabel") or stanza["track"],
                    "assemblyNames": [genome],
                    "category": ["HPRC", stanza.get("shortLabel", "Annotations")],
                    "adapter": {
                        "type": "BigBedAdapter",
                        "bigBedLocation": {
                            "uri": uri,
                            "locationType": "UriLocation",
                        },
                    },
                }
            )

    config = {
        "assemblies": list(assemblies.values()),
        "tracks": tracks,
        "defaultSession": {
            "name": "HPRC",
            "views": [],
        },
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {out} with {len(assemblies)} assemblies and {len(tracks)} tracks")


if __name__ == "__main__":
    main()
