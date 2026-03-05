process AGGREGATE_TRACKHUB {
    tag "$aggregate_name"
    label 'process_low'

    container "${ params.container_python ?: 'python:3.11-slim' }"

    publishDir "${params.outdir}/trackhubs", mode: 'copy'

    input:
    path manifest_json
    val aggregate_name
    val outdir
    val email
    val short_label
    val long_label

    output:
    path "${aggregate_name}/**", emit: hub

    script:
    """
    python - <<'PY'
import os, json, sys

with open('${manifest_json}', 'r') as fh:
    hubs = json.load(fh)

agg_name = "${aggregate_name}"
outdir = "${outdir}"
email = "${email}"
short = "${short_label}" if "${short_label}" else agg_name[:17]
long = "${long_label}" if "${long_label}" else f"{agg_name} multi-assembly hub"

root = os.path.join(outdir, 'trackhubs')
agg_dir = os.path.join(agg_name)
os.makedirs(agg_dir, exist_ok=True)

# hub.txt
with open(os.path.join(agg_dir, 'hub.txt'), 'w') as fh:
    fh.write(f"hub {agg_name}\n")
    fh.write(f"shortLabel {short}\n")
    fh.write(f"longLabel {long}\n")
    fh.write(f"email {email}\n")

# genomes.txt
with open(os.path.join(agg_dir, 'genomes.txt'), 'w') as fh:
    for entry in hubs:
        genome = entry['genome']
        trackdb_abs = entry['trackdb']
        target = os.path.join(root, agg_name)
        rel = os.path.relpath(trackdb_abs, start=target)
        fh.write(f"genome {genome}\n")
        fh.write(f"trackDb {rel}\n\n")

PY
    """

    stub:
    """
    mkdir -p ${aggregate_name}
    printf 'stub' > ${aggregate_name}/hub.txt
    printf 'stub' > versions.yml
    """
}
