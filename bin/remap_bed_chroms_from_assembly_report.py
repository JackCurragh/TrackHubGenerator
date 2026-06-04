#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from typing import Dict, Iterable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Remap BED chromosome names using aliases from an NCBI assembly_report.txt. "
            "Rows that cannot be mapped are left unchanged so downstream validation can fail loudly."
        )
    )
    parser.add_argument("--report", required=True, help="NCBI assembly_report.txt")
    parser.add_argument("--target-column", required=True, help="Target assembly report column, e.g. genbank-accn")
    parser.add_argument("--input", required=True, help="Input BED")
    parser.add_argument("--output", required=True, help="Output BED")
    return parser.parse_args()


def norm_header(name: str) -> str:
    return name.strip().lower().replace(" ", "-").replace("_", "-")


def read_headers(report: str) -> list[str]:
    headers: list[str] = []
    with open(report, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("#"):
                continue
            if "Sequence-Name" in line and ("Sequence-Length" in line or "sequence length" in line.lower()):
                text = line.lstrip("#").strip()
                if "\t" not in text:
                    text = "\t".join(text.split())
                headers = [norm_header(value) for value in text.split("\t")]
    return headers


def non_empty_aliases(parts: list[str], indexes: Iterable[int]) -> list[str]:
    values: list[str] = []
    for index in indexes:
        if index >= len(parts):
            continue
        value = parts[index].strip()
        if value and value.lower() != "na":
            values.append(value)
    return values


def build_alias_map(report: str, target_column: str) -> Dict[str, str]:
    headers = read_headers(report)
    if not headers:
        raise SystemExit(f"Could not parse assembly report header from {report}")

    idx = {name: i for i, name in enumerate(headers)}
    target_column = norm_header(target_column)
    if target_column not in idx:
        raise SystemExit(f"Target column {target_column!r} is not present in {report}")

    alias_columns = [
        "sequence-name",
        "assigned-molecule",
        "genbank-accn",
        "refseq-accn",
        "ucsc-style-name",
    ]
    alias_indexes = [idx[name] for name in alias_columns if name in idx]
    target_index = idx[target_column]

    aliases: Dict[str, str] = {}
    with open(report, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if target_index >= len(parts):
                continue
            target = parts[target_index].strip()
            if not target or target.lower() == "na":
                continue
            for alias in non_empty_aliases(parts, alias_indexes):
                aliases[alias] = target
                if alias == "MT":
                    aliases["M"] = target
                elif alias == "M":
                    aliases["MT"] = target
    return aliases


def main() -> None:
    args = parse_args()
    aliases = build_alias_map(args.report, args.target_column)
    changed = 0
    total = 0

    with open(args.input, "rt", encoding="utf-8", errors="ignore") as src, open(args.output, "wt") as dst:
        for line in src:
            if not line.strip() or line.startswith("#"):
                dst.write(line)
                continue
            total += 1
            parts = line.rstrip("\n").split("\t")
            mapped = aliases.get(parts[0], parts[0])
            if mapped != parts[0]:
                changed += 1
                parts[0] = mapped
            dst.write("\t".join(parts) + "\n")

    print(f"Remapped {changed}/{total} BED rows using {len(aliases)} assembly-report aliases", file=sys.stderr)


if __name__ == "__main__":
    main()
