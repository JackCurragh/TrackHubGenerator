#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import tempfile
from collections import defaultdict
from collections import OrderedDict
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a static per-gene lookup from the GFF inputs used to create HPRC track hubs. "
            "Coordinates are normalized to GenBank accessions with the NCBI assembly reports."
        )
    )
    parser.add_argument("--input-csv", required=True, help="Trackhub input CSV with Run,Path,Genome,AssemblyReport,Hub")
    parser.add_argument("--trackhubs-root", required=True, help="Generated trackhubs directory")
    parser.add_argument("--output-dir", required=True, help="Output lookup directory")
    parser.add_argument(
        "--base-url",
        default="",
        help="Optional public URL corresponding to --trackhubs-root; used for hub/BigBed links",
    )
    parser.add_argument(
        "--tmp-dir",
        default="",
        help="Directory for temporary JSONL shards. Defaults to a temporary directory beside --output-dir.",
    )
    parser.add_argument(
        "--shards",
        type=int,
        default=1024,
        help="Number of temporary hash shards used while streaming GFF entries.",
    )
    parser.add_argument(
        "--keep-tmp",
        action="store_true",
        help="Keep temporary shards after a successful run for debugging.",
    )
    return parser.parse_args()


def parse_attrs(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in text.rstrip().split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip('"')
        else:
            continue
        attrs[key.strip()] = value.strip()
    return attrs


def pick_attr(attrs: dict[str, str], keys: list[str]) -> str:
    for key in keys:
        value = attrs.get(key)
        if value and value != "NA":
            return value
    return ""


def sample_haplotype_from_hub(hub: str) -> tuple[str, str]:
    name = hub.removeprefix("HPRC_")
    match = re.match(r"^(?P<sample>.+)_(?P<hap>mat|pat|hap1|hap2)$", name)
    if not match:
        return name, ""
    hap = match.group("hap")
    hap = {"hap1": "pat", "hap2": "mat"}.get(hap, hap)
    return match.group("sample"), hap


def source_from_run(run: str) -> str:
    if run.endswith("_CAT"):
        return "CAT"
    if run.endswith("_ENSEMBL"):
        return "ENSEMBL"
    return run.rsplit("_", 1)[-1]


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("rt", encoding="utf-8", errors="ignore")


def normalise_header(name: str) -> str:
    return name.strip().lower().replace(" ", "-").replace("_", "-")


def assembly_report_aliases(report: Path) -> dict[str, str]:
    headers: list[str] = []
    with report.open(errors="ignore") as handle:
        for line in handle:
            if line.startswith("#") and "Sequence-Name" in line and "Sequence-Length" in line:
                text = line.lstrip("#").strip()
                if "\t" not in text:
                    text = "\t".join(text.split())
                headers = [normalise_header(x) for x in text.split("\t")]
                break

    aliases: dict[str, str] = {}
    if not headers:
        return aliases

    idx = {name: i for i, name in enumerate(headers)}
    target_i = idx.get("genbank-accn")
    if target_i is None:
        return aliases

    alias_columns = [
        "sequence-name",
        "assigned-molecule",
        "genbank-accn",
        "refseq-accn",
        "ucsc-style-name",
    ]
    alias_indexes = [idx[name] for name in alias_columns if name in idx]

    with report.open(errors="ignore") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if target_i >= len(fields):
                continue
            target = fields[target_i]
            if not target or target == "na":
                continue
            for alias_i in alias_indexes:
                if alias_i < len(fields):
                    alias = fields[alias_i]
                    if alias and alias != "na":
                        aliases[alias] = target
            assigned_i = idx.get("assigned-molecule")
            if assigned_i is not None and assigned_i < len(fields):
                assigned = fields[assigned_i]
                if assigned == "MT":
                    aliases["M"] = target
                elif assigned == "M":
                    aliases["MT"] = target
    return aliases


def safe_gene_filename(key: str) -> tuple[str, str]:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", key).strip("_") or "unknown"
    prefix = safe[0].upper() if safe[0].isalnum() else "_"
    return prefix, safe


def public_url(base_url: str, root: Path, path: Path) -> str:
    if not base_url:
        return ""
    rel = path.resolve().relative_to(root.resolve()).as_posix()
    return f"{base_url.rstrip('/')}/{rel}"


def gene_key(symbol: str, gene_id: str) -> str:
    return symbol or gene_id or "unknown"


def shard_for_key(key: str, n_shards: int) -> int:
    digest = hashlib.blake2s(key.encode("utf-8"), digest_size=4).hexdigest()
    return int(digest, 16) % n_shards


def write_jsonl(handle, key: str, entry: dict[str, Any]) -> None:
    handle.write(json.dumps({"key": key, "entry": entry}, separators=(",", ":")))
    handle.write("\n")


class ShardWriters:
    def __init__(self, shard_paths: list[Path], max_open: int = 64) -> None:
        self.shard_paths = shard_paths
        self.max_open = max_open
        self.handles: OrderedDict[int, Any] = OrderedDict()

    def write(self, shard: int, key: str, entry: dict[str, Any]) -> None:
        handle = self.handles.get(shard)
        if handle is None:
            handle = self.shard_paths[shard].open("a", encoding="utf-8")
            self.handles[shard] = handle
            if len(self.handles) > self.max_open:
                _, old_handle = self.handles.popitem(last=False)
                old_handle.close()
        else:
            self.handles.move_to_end(shard)
        write_jsonl(handle, key, entry)

    def close(self) -> None:
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()

    def __enter__(self) -> "ShardWriters":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


def tmp_root_for(args: argparse.Namespace, output_dir: Path) -> Path:
    if args.tmp_dir:
        path = Path(args.tmp_dir)
        path.mkdir(parents=True, exist_ok=True)
        return Path(tempfile.mkdtemp(prefix="gene_lookup_", dir=path))
    output_dir.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix=".gene_lookup_", dir=output_dir))


def main() -> None:
    args = parse_args()
    if args.shards < 1:
        raise SystemExit("--shards must be >= 1")

    input_csv = Path(args.input_csv)
    trackhubs_root = Path(args.trackhubs_root)
    output_dir = Path(args.output_dir)
    genes_dir = output_dir / "genes"
    output_dir.mkdir(parents=True, exist_ok=True)
    genes_dir.mkdir(parents=True, exist_ok=True)
    tmp_root = tmp_root_for(args, output_dir)
    shard_paths = [tmp_root / f"shard_{i:04d}.jsonl" for i in range(args.shards)]

    assemblies: dict[str, dict[str, Any]] = {}
    report_cache: dict[Path, dict[str, str]] = {}
    stats = {
        "rows": 0,
        "gene_entries": 0,
        "missing_bigbeds": 0,
        "unmapped_seqids": 0,
    }

    try:
        with input_csv.open(newline="") as handle, ShardWriters(shard_paths) as shard_writers:
            for row in csv.DictReader(handle):
                stats["rows"] += 1
                run = row["Run"]
                genome = row["Genome"]
                hub = row["Hub"]
                source = source_from_run(run)
                sample, haplotype = sample_haplotype_from_hub(hub)
                gff = Path(row["Path"])
                report = Path(row["AssemblyReport"])
                bigbed = trackhubs_root / hub / genome / f"{hub}_{run}.bigBed"
                trackdb = trackhubs_root / hub / genome / "trackDb.txt"
                hub_file = trackhubs_root / hub / "hub.txt"

                assembly = assemblies.setdefault(
                    genome,
                    {
                        "assembly_accession": genome,
                        "sample": sample,
                        "haplotype": haplotype,
                        "hub": hub,
                        "hub_path": str(hub_file),
                        "hub_url": public_url(args.base_url, trackhubs_root, hub_file) if hub_file.exists() else "",
                        "trackdb_path": str(trackdb),
                        "tracks": {},
                    },
                )
                assembly["tracks"][source] = {
                    "run": run,
                    "source_gff": str(gff),
                    "bigbed_path": str(bigbed),
                    "bigbed_url": public_url(args.base_url, trackhubs_root, bigbed) if bigbed.exists() else "",
                }

                if not bigbed.exists():
                    stats["missing_bigbeds"] += 1

                aliases = report_cache.setdefault(report, assembly_report_aliases(report))

                with open_text(gff) as gff_handle:
                    for line in gff_handle:
                        if line.startswith("#"):
                            continue
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) < 9 or fields[2] != "gene":
                            continue
                        seqid = fields[0]
                        genbank_seqid = aliases.get(seqid, seqid)
                        if genbank_seqid == seqid and aliases and seqid not in aliases.values():
                            stats["unmapped_seqids"] += 1

                        attrs = parse_attrs(fields[8])
                        gene_id = pick_attr(attrs, ["gene_id", "ID", "Dbxref"])
                        symbol = pick_attr(attrs, ["gene_name", "Name", "gene", "standard_name"])
                        name = pick_attr(attrs, ["description", "Note", "product"])
                        biotype = pick_attr(attrs, ["gene_biotype", "biotype", "gene_type"])
                        key = gene_key(symbol, gene_id)

                        entry = {
                            "sample": sample,
                            "haplotype": haplotype,
                            "assembly_accession": genome,
                            "hub": hub,
                            "source": source,
                            "gene_id": gene_id,
                            "symbol": symbol,
                            "name": name,
                            "biotype": biotype,
                            "region_name": seqid,
                            "genbank_seqid": genbank_seqid,
                            "start": int(fields[3]),
                            "end": int(fields[4]),
                            "strand": fields[6],
                            "trackhub_path": str(hub_file),
                            "trackhub_url": public_url(args.base_url, trackhubs_root, hub_file) if hub_file.exists() else "",
                            "bigbed_path": str(bigbed),
                            "bigbed_url": public_url(args.base_url, trackhubs_root, bigbed) if bigbed.exists() else "",
                        }
                        shard_writers.write(shard_for_key(key, args.shards), key, entry)
                        stats["gene_entries"] += 1
                if stats["rows"] % 25 == 0:
                    print(f"Processed {stats['rows']} input rows and {stats['gene_entries']} gene entries", flush=True)

        gene_index = []
        for shard_path in shard_paths:
            if not shard_path.exists():
                continue
            genes: dict[str, list[dict[str, Any]]] = defaultdict(list)
            with shard_path.open(encoding="utf-8") as shard_handle:
                for line in shard_handle:
                    item = json.loads(line)
                    genes[item["key"]].append(item["entry"])

            for key, entries in sorted(genes.items(), key=lambda item: item[0].lower()):
                prefix, safe = safe_gene_filename(key)
                rel_path = Path("genes") / prefix / f"{safe}.json.gz"
                out_path = output_dir / rel_path
                out_path.parent.mkdir(parents=True, exist_ok=True)
                symbols = sorted({e["symbol"] for e in entries if e["symbol"]})
                gene_ids = sorted({e["gene_id"] for e in entries if e["gene_id"]})
                sources = sorted({e["source"] for e in entries})
                samples = sorted({e["sample"] for e in entries})
                payload = {
                    "key": key,
                    "symbols": symbols,
                    "gene_ids": gene_ids,
                    "entries": sorted(
                        entries,
                        key=lambda e: (e["sample"], e["haplotype"], e["assembly_accession"], e["source"], e["start"]),
                    ),
                }
                with gzip.open(out_path, "wt", encoding="utf-8") as out_handle:
                    json.dump(payload, out_handle, separators=(",", ":"))
                gene_index.append(
                    {
                        "key": key,
                        "symbols": symbols,
                        "gene_ids": gene_ids[:20],
                        "sources": sources,
                        "n_entries": len(entries),
                        "n_samples": len(samples),
                        "path": rel_path.as_posix(),
                    }
                )

        gene_index.sort(key=lambda item: item["key"].lower())

        with (output_dir / "assemblies.json").open("w", encoding="utf-8") as handle:
            json.dump(sorted(assemblies.values(), key=lambda x: (x["sample"], x["haplotype"])), handle, indent=2)

        with gzip.open(output_dir / "gene_index.json.gz", "wt", encoding="utf-8") as handle:
            json.dump(gene_index, handle, separators=(",", ":"))

        manifest = {
            "assemblies": len(assemblies),
            "genes": len(gene_index),
            "tmp_dir": str(tmp_root) if args.keep_tmp else "",
            **stats,
            "files": {
                "assemblies": "assemblies.json",
                "gene_index": "gene_index.json.gz",
                "genes": "genes/<prefix>/<gene>.json.gz",
            },
        }
        with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2)

        print(json.dumps(manifest, indent=2))
    finally:
        if not args.keep_tmp:
            shutil.rmtree(tmp_root, ignore_errors=True)


if __name__ == "__main__":
    main()
