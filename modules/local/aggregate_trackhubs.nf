process AGGREGATE_TRACKHUB {
    tag "$aggregate_name"
    label 'process_low'

    container "${ params.container_python ?: 'python:3.11-slim' }"

    publishDir "${params.outdir}/trackhubs", mode: 'copy'

    input:
    path writer_script
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
    chmod +x ${writer_script}
    ./${writer_script.getName()} \
        --manifest '${manifest_json}' \
        --aggregate-name "${aggregate_name}" \
        --outdir "${outdir}" \
        --email "${email}" \
        --short-label "${short_label}" \
        --long-label "${long_label}"
    """

    stub:
    """
    mkdir -p ${aggregate_name}
    printf 'stub' > ${aggregate_name}/hub.txt
    printf 'stub' > versions.yml
    """
}
