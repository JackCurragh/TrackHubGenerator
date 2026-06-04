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
import time
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed
from collections import defaultdict
from collections import OrderedDict
from pathlib import Path
from typing import Any


def add_shared_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--base-url",
        default="",
        help="Optional public URL corresponding to --trackhubs-root; used for hub/BigBed links",
    )
    parser.add_argument(
        "--shards",
        type=int,
        default=1024,
        help="Number of temporary hash shards used while streaming GFF entries.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=10,
        help="Print progress every N input CSV rows or N shard files. Use 0 to disable progress.",
    )
    parser.add_argument(
        "--include-unnamed",
        action="store_true",
        help=(
            "Include gene rows with no usable symbol. By default these are skipped because CAT can contain "
            "millions of assembly-local IDs that make the browser lookup large and hard to use."
        ),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build a static per-gene lookup from the GFF inputs used to create HPRC track hubs. "
            "Coordinates are normalized to GenBank accessions with the NCBI assembly reports."
        )
    )
    parser.add_argument(
        "command",
        nargs="?",
        choices=["build", "scan-row", "merge"],
        default="build",
        help="Mode to run. Omit for the standard single-command build.",
    )
    parser.add_argument("--input-csv", required=True, help="Trackhub input CSV with Run,Path,Genome,AssemblyReport,Hub")
    parser.add_argument("--trackhubs-root", required=True, help="Generated trackhubs directory")
    parser.add_argument("--output-dir", required=True, help="Output lookup directory")
    parser.add_argument(
        "--tmp-dir",
        default="",
        help="Directory for temporary JSONL shards. Defaults to a temporary directory beside --output-dir.",
    )
    parser.add_argument(
        "--keep-tmp",
        action="store_true",
        help="Keep temporary shards after a successful run for debugging.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel GFF input rows to scan. Use 1 for serial execution.",
    )
    parser.add_argument(
        "--row",
        type=int,
        help="1-based input CSV row to scan in scan-row mode. Usually SLURM_ARRAY_TASK_ID.",
    )
    parser.add_argument(
        "--parts-dir",
        default="",
        help="Directory containing scan-row shard outputs for merge mode.",
    )
    add_shared_arguments(parser)
    return parser


def parse_args() -> argparse.Namespace:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "scan-row" and not args.row:
        parser.error("scan-row requires --row")
    if args.command == "merge" and not args.parts_dir:
        parser.error("merge requires --parts-dir")
    return args


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


def is_useful_lookup_gene(source: str, symbol: str, gene_id: str, include_unnamed: bool) -> bool:
    if include_unnamed:
        return bool(symbol or gene_id)
    if symbol:
        return True
    return source == "ENSEMBL" and gene_id.startswith("ENS")


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


def count_csv_rows(path: Path) -> int:
    with path.open(newline="") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def format_elapsed(seconds: float) -> str:
    seconds = int(seconds)
    hours, rem = divmod(seconds, 3600)
    minutes, seconds = divmod(rem, 60)
    if hours:
        return f"{hours}h{minutes:02d}m{seconds:02d}s"
    if minutes:
        return f"{minutes}m{seconds:02d}s"
    return f"{seconds}s"


def print_scan_progress(
    done: int,
    total: int,
    genes: int,
    skipped: int,
    start_time: float,
    label: str = "",
) -> None:
    elapsed = time.monotonic() - start_time
    pct = 100 * done / total if total else 100
    rate = 60 * done / elapsed if elapsed > 0 else 0
    suffix = f" {label}" if label else ""
    print(
        f"[scan] {done}/{total} rows ({pct:.1f}%) | indexed={genes:,} | skipped_unnamed={skipped:,} | "
        f"elapsed={format_elapsed(elapsed)} | rate={rate:.2f} rows/min{suffix}",
        flush=True,
    )


def print_write_progress(done: int, total: int, genes: int, start_time: float) -> None:
    elapsed = time.monotonic() - start_time
    pct = 100 * done / total if total else 100
    rate = done / elapsed if elapsed > 0 else 0
    print(
        f"[write] {done}/{total} shards ({pct:.1f}%) | genes={genes:,} | "
        f"elapsed={format_elapsed(elapsed)} | rate={rate:.2f} shards/sec",
        flush=True,
    )


def empty_scan_stats() -> dict[str, int]:
    return {
        "rows": 0,
        "gene_features_seen": 0,
        "gene_entries": 0,
        "skipped_unnamed": 0,
        "missing_bigbeds": 0,
        "unmapped_seqids": 0,
    }


def add_scan_stats(target: dict[str, int], source: dict[str, int]) -> None:
    for key, value in source.items():
        target[key] = target.get(key, 0) + value


def scan_row_to_shards(
    row_i: int,
    row: dict[str, str],
    trackhubs_root_text: str,
    tmp_root_text: str,
    n_shards: int,
    base_url: str,
    include_unnamed: bool,
) -> dict[str, Any]:
    trackhubs_root = Path(trackhubs_root_text)
    tmp_root = Path(tmp_root_text)
    run = row["Run"]
    genome = row["Genome"]
    hub = row["Hub"]
    source = source_from_run(run)
    sample, haplotype = sample_haplotype_from_hub(hub)
    gff = Path(row["Path"])
    report = Path(row["AssemblyReport"])
    bigbed = trackhubs_root / hub / genome / f"{hub}_{run}.bigBed"
    hub_file = trackhubs_root / hub / "hub.txt"
    worker_dir = tmp_root / f"worker_{row_i:06d}"
    worker_dir.mkdir(parents=True, exist_ok=True)
    shard_paths = [worker_dir / f"shard_{i:04d}.jsonl" for i in range(n_shards)]
    aliases = assembly_report_aliases(report)
    stats = empty_scan_stats()
    stats["rows"] = 1
    if not bigbed.exists():
        stats["missing_bigbeds"] += 1

    with ShardWriters(shard_paths) as shard_writers, open_text(gff) as gff_handle:
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
            stats["gene_features_seen"] += 1
            if not is_useful_lookup_gene(source, symbol, gene_id, include_unnamed):
                stats["skipped_unnamed"] += 1
                continue
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
                "trackhub_url": public_url(base_url, trackhubs_root, hub_file) if hub_file.exists() else "",
                "bigbed_path": str(bigbed),
                "bigbed_url": public_url(base_url, trackhubs_root, bigbed) if bigbed.exists() else "",
            }
            shard_writers.write(shard_for_key(key, n_shards), key, entry)
            stats["gene_entries"] += 1

    return {"row_i": row_i, "run": run, "stats": stats}


def add_assembly_metadata(
    assemblies: dict[str, dict[str, Any]],
    row: dict[str, str],
    trackhubs_root: Path,
    base_url: str,
) -> None:
    run = row["Run"]
    genome = row["Genome"]
    hub = row["Hub"]
    source = source_from_run(run)
    sample, haplotype = sample_haplotype_from_hub(hub)
    gff = Path(row["Path"])
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
            "hub_url": public_url(base_url, trackhubs_root, hub_file) if hub_file.exists() else "",
            "trackdb_path": str(trackdb),
            "tracks": {},
        },
    )
    assembly["tracks"][source] = {
        "run": run,
        "source_gff": str(gff),
        "bigbed_path": str(bigbed),
        "bigbed_url": public_url(base_url, trackhubs_root, bigbed) if bigbed.exists() else "",
    }


def read_input_rows(input_csv: Path) -> list[dict[str, str]]:
    with input_csv.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_scan_result(parts_dir: Path, result: dict[str, Any]) -> None:
    stats_dir = parts_dir / "stats"
    stats_dir.mkdir(parents=True, exist_ok=True)
    row_i = int(result["row_i"])
    with (stats_dir / f"row_{row_i:06d}.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)


def read_scan_stats(parts_dir: Path, expected_rows: int) -> dict[str, int]:
    stats = empty_scan_stats()
    stats_dir = parts_dir / "stats"
    seen_rows = set()
    for path in sorted(stats_dir.glob("row_*.json")):
        with path.open(encoding="utf-8") as handle:
            result = json.load(handle)
        seen_rows.add(int(result["row_i"]))
        add_scan_stats(stats, result["stats"])

    missing = sorted(set(range(1, expected_rows + 1)) - seen_rows)
    if missing:
        preview = ",".join(str(x) for x in missing[:20])
        suffix = "..." if len(missing) > 20 else ""
        raise SystemExit(f"Missing scan-row outputs for {len(missing)} row(s): {preview}{suffix}")
    return stats


def scan_rows(
    rows: list[dict[str, str]],
    trackhubs_root: Path,
    parts_dir: Path,
    base_url: str,
    shards: int,
    include_unnamed: bool,
    jobs: int,
    progress_every: int,
) -> dict[str, int]:
    stats = empty_scan_stats()
    start_time = time.monotonic()
    print(
        f"Scanning {len(rows)} input rows with {jobs} job(s) into {shards} temporary shards at {parts_dir}",
        flush=True,
    )

    if jobs == 1:
        for row_i, row in enumerate(rows, start=1):
            result = scan_row_to_shards(row_i, row, str(trackhubs_root), str(parts_dir), shards, base_url, include_unnamed)
            write_scan_result(parts_dir, result)
            add_scan_stats(stats, result["stats"])
            if progress_every and stats["rows"] % progress_every == 0:
                print_scan_progress(
                    stats["rows"],
                    len(rows),
                    stats["gene_entries"],
                    stats["skipped_unnamed"],
                    start_time,
                    result["run"],
                )
    else:
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            futures = [
                executor.submit(
                    scan_row_to_shards,
                    row_i,
                    row,
                    str(trackhubs_root),
                    str(parts_dir),
                    shards,
                    base_url,
                    include_unnamed,
                )
                for row_i, row in enumerate(rows, start=1)
            ]
            for future in as_completed(futures):
                result = future.result()
                write_scan_result(parts_dir, result)
                add_scan_stats(stats, result["stats"])
                if progress_every and stats["rows"] % progress_every == 0:
                    print_scan_progress(
                        stats["rows"],
                        len(rows),
                        stats["gene_entries"],
                        stats["skipped_unnamed"],
                        start_time,
                        result["run"],
                    )

    if not progress_every or stats["rows"] % progress_every != 0:
        print_scan_progress(
            stats["rows"],
            len(rows),
            stats["gene_entries"],
            stats["skipped_unnamed"],
            start_time,
            "complete",
        )
    return stats


def build_assemblies(rows: list[dict[str, str]], trackhubs_root: Path, base_url: str) -> dict[str, dict[str, Any]]:
    assemblies: dict[str, dict[str, Any]] = {}
    for row in rows:
        add_assembly_metadata(assemblies, row, trackhubs_root, base_url)
    return assemblies


def collect_shard_groups(parts_dir: Path) -> dict[int, list[Path]]:
    shard_groups: dict[int, list[Path]] = defaultdict(list)
    for worker_dir in parts_dir.glob("worker_*"):
        for shard_path in worker_dir.glob("shard_*.jsonl"):
            shard_id = int(shard_path.stem.split("_")[1])
            shard_groups[shard_id].append(shard_path)
    return shard_groups


def merge_parts(
    rows: list[dict[str, str]],
    trackhubs_root: Path,
    parts_dir: Path,
    output_dir: Path,
    base_url: str,
    progress_every: int,
    tmp_label: str = "",
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    genes_dir = output_dir / "genes"
    genes_dir.mkdir(parents=True, exist_ok=True)
    assemblies = build_assemblies(rows, trackhubs_root, base_url)
    stats = read_scan_stats(parts_dir, len(rows))

    gene_index = []
    shard_groups = collect_shard_groups(parts_dir)
    existing_shards = sorted(shard_groups)
    write_start = time.monotonic()
    print(f"Writing per-gene JSON from {len(existing_shards)} non-empty shards", flush=True)
    for shard_i, shard_id in enumerate(existing_shards, start=1):
        genes: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for shard_path in shard_groups[shard_id]:
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
        if progress_every and shard_i % progress_every == 0:
            print_write_progress(shard_i, len(existing_shards), len(gene_index), write_start)
    if not progress_every or len(existing_shards) % progress_every != 0:
        print_write_progress(len(existing_shards), len(existing_shards), len(gene_index), write_start)

    gene_index.sort(key=lambda item: item["key"].lower())

    with (output_dir / "assemblies.json").open("w", encoding="utf-8") as handle:
        json.dump(sorted(assemblies.values(), key=lambda x: (x["sample"], x["haplotype"])), handle, indent=2)

    with gzip.open(output_dir / "gene_index.json.gz", "wt", encoding="utf-8") as handle:
        json.dump(gene_index, handle, separators=(",", ":"))

    manifest = {
        "assemblies": len(assemblies),
        "genes": len(gene_index),
        "tmp_dir": tmp_label,
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
    return manifest


def run_build(args: argparse.Namespace) -> None:
    input_csv = Path(args.input_csv)
    trackhubs_root = Path(args.trackhubs_root)
    output_dir = Path(args.output_dir)
    rows = read_input_rows(input_csv)
    tmp_root = tmp_root_for(args, output_dir)
    try:
        scan_rows(
            rows,
            trackhubs_root,
            tmp_root,
            args.base_url,
            args.shards,
            args.include_unnamed,
            args.jobs,
            args.progress_every,
        )
        merge_parts(
            rows,
            trackhubs_root,
            tmp_root,
            output_dir,
            args.base_url,
            args.progress_every,
            str(tmp_root) if args.keep_tmp else "",
        )
    finally:
        if not args.keep_tmp:
            shutil.rmtree(tmp_root, ignore_errors=True)


def run_scan_row(args: argparse.Namespace) -> None:
    rows = read_input_rows(Path(args.input_csv))
    if args.row < 1 or args.row > len(rows):
        raise SystemExit(f"--row must be between 1 and {len(rows)}")
    parts_dir = Path(args.output_dir)
    parts_dir.mkdir(parents=True, exist_ok=True)
    result = scan_row_to_shards(
        args.row,
        rows[args.row - 1],
        args.trackhubs_root,
        str(parts_dir),
        args.shards,
        args.base_url,
        args.include_unnamed,
    )
    write_scan_result(parts_dir, result)
    print(json.dumps(result, indent=2), flush=True)


def run_merge(args: argparse.Namespace) -> None:
    rows = read_input_rows(Path(args.input_csv))
    merge_parts(
        rows,
        Path(args.trackhubs_root),
        Path(args.parts_dir),
        Path(args.output_dir),
        args.base_url,
        args.progress_every,
        "",
    )


def main() -> None:
    args = parse_args()
    if args.shards < 1:
        raise SystemExit("--shards must be >= 1")
    if args.jobs < 1:
        raise SystemExit("--jobs must be >= 1")

    if args.command == "scan-row":
        run_scan_row(args)
    elif args.command == "merge":
        run_merge(args)
    else:
        run_build(args)


if __name__ == "__main__":
    main()
