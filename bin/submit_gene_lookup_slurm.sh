#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Submit HPRC gene lookup as a Slurm scan array plus dependent merge job.

Required:
  --input-csv PATH
  --trackhubs-root PATH
  --lookup-dir PATH
  --parts-dir PATH
  --base-url URL

Optional:
  --shards N              default: 4096
  --array-concurrency N   default: 50
  --scan-mem MEM          default: 4G
  --scan-time TIME        default: 02:00:00
  --merge-mem MEM         default: 16G
  --merge-time TIME       default: 04:00:00
  --progress-every N      default: 25
  --include-unnamed       include unnamed/local gene IDs

Example:
  bin/submit_gene_lookup_slurm.sh \
    --input-csv "$OUT/hprc_trackhub_input.csv" \
    --trackhubs-root "$OUT/trackhubs" \
    --lookup-dir "$OUT/trackhubs/lookup" \
    --parts-dir "$OUT/gene_lookup_parts" \
    --base-url "https://ftp.ebi.ac.uk/pub/databases/ensembl/hprc/release_2/R2_trackhubs"
EOF
}

INPUT_CSV=""
TRACKHUBS_ROOT=""
LOOKUP_DIR=""
PARTS_DIR=""
BASE_URL=""
SHARDS=4096
ARRAY_CONCURRENCY=50
SCAN_MEM=4G
SCAN_TIME=02:00:00
MERGE_MEM=16G
MERGE_TIME=04:00:00
PROGRESS_EVERY=25
INCLUDE_UNNAMED=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-csv) INPUT_CSV="$2"; shift 2 ;;
    --trackhubs-root) TRACKHUBS_ROOT="$2"; shift 2 ;;
    --lookup-dir) LOOKUP_DIR="$2"; shift 2 ;;
    --parts-dir) PARTS_DIR="$2"; shift 2 ;;
    --base-url) BASE_URL="$2"; shift 2 ;;
    --shards) SHARDS="$2"; shift 2 ;;
    --array-concurrency) ARRAY_CONCURRENCY="$2"; shift 2 ;;
    --scan-mem) SCAN_MEM="$2"; shift 2 ;;
    --scan-time) SCAN_TIME="$2"; shift 2 ;;
    --merge-mem) MERGE_MEM="$2"; shift 2 ;;
    --merge-time) MERGE_TIME="$2"; shift 2 ;;
    --progress-every) PROGRESS_EVERY="$2"; shift 2 ;;
    --include-unnamed) INCLUDE_UNNAMED=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for name in INPUT_CSV TRACKHUBS_ROOT LOOKUP_DIR PARTS_DIR BASE_URL; do
  if [[ -z "${!name}" ]]; then
    echo "Missing required argument for $name" >&2
    usage >&2
    exit 2
  fi
done

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
LOOKUP_SCRIPT="$SCRIPT_DIR/build_gene_lookup.py"
N_ROWS=$(($(wc -l < "$INPUT_CSV") - 1))
if [[ "$N_ROWS" -lt 1 ]]; then
  echo "No input rows found in $INPUT_CSV" >&2
  exit 1
fi

mkdir -p "$PARTS_DIR/logs" "$LOOKUP_DIR"
SCAN_SCRIPT="$PARTS_DIR/scan_gene_lookup.sbatch"
MERGE_SCRIPT="$PARTS_DIR/merge_gene_lookup.sbatch"
EXTRA_ARGS=()
if [[ "$INCLUDE_UNNAMED" -eq 1 ]]; then
  EXTRA_ARGS+=(--include-unnamed)
fi

cat > "$SCAN_SCRIPT" <<EOF
#!/usr/bin/env bash
#SBATCH -J hprc_gene_lookup_scan
#SBATCH --mem=$SCAN_MEM
#SBATCH --time=$SCAN_TIME
#SBATCH -o $PARTS_DIR/logs/scan_%A_%a.out
#SBATCH -e $PARTS_DIR/logs/scan_%A_%a.err

set -euo pipefail

"$LOOKUP_SCRIPT" scan-row \\
  --input-csv "$INPUT_CSV" \\
  --trackhubs-root "$TRACKHUBS_ROOT" \\
  --output-dir "$PARTS_DIR" \\
  --base-url "$BASE_URL" \\
  --shards "$SHARDS" \\
  --row "\$SLURM_ARRAY_TASK_ID" ${EXTRA_ARGS[*]+"${EXTRA_ARGS[*]}"}
EOF

cat > "$MERGE_SCRIPT" <<EOF
#!/usr/bin/env bash
#SBATCH -J hprc_gene_lookup_merge
#SBATCH --mem=$MERGE_MEM
#SBATCH --time=$MERGE_TIME
#SBATCH -o $PARTS_DIR/logs/merge_%j.out
#SBATCH -e $PARTS_DIR/logs/merge_%j.err

set -euo pipefail

"$LOOKUP_SCRIPT" merge \\
  --input-csv "$INPUT_CSV" \\
  --trackhubs-root "$TRACKHUBS_ROOT" \\
  --parts-dir "$PARTS_DIR" \\
  --output-dir "$LOOKUP_DIR" \\
  --base-url "$BASE_URL" \\
  --shards "$SHARDS" \\
  --progress-every "$PROGRESS_EVERY" ${EXTRA_ARGS[*]+"${EXTRA_ARGS[*]}"}
EOF

scan_job=$(sbatch --parsable --array=1-"$N_ROWS"%"$ARRAY_CONCURRENCY" "$SCAN_SCRIPT")
merge_job=$(sbatch --parsable --dependency=afterok:"$scan_job" "$MERGE_SCRIPT")

cat <<EOF
Submitted gene lookup scan array: $scan_job
Submitted dependent merge job:    $merge_job

Input rows: $N_ROWS
Parts dir:  $PARTS_DIR
Lookup dir: $LOOKUP_DIR
Logs:       $PARTS_DIR/logs
EOF
