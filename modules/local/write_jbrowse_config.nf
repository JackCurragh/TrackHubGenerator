process WRITE_JBROWSE_CONFIG {
    tag "jbrowse_config"
    label 'process_low'

    container "${ params.container_python ?: 'python:3.11-slim' }"

    publishDir "${params.outdir}/jbrowse", mode: 'copy'

    input:
    path writer_script
    path trackhubs_root
    val base_url
    val twobit_url_template
    val chrom_sizes_url_template

    output:
    path "config.json", emit: config

    script:
    def chromSizesArg = chrom_sizes_url_template ? "--chrom-sizes-url-template '${chrom_sizes_url_template}'" : ''
    """
    chmod +x ${writer_script}
    ./${writer_script.getName()} \
        --trackhubs-root '${trackhubs_root}' \
        --base-url '${base_url}' \
        --twobit-url-template '${twobit_url_template}' \
        ${chromSizesArg} \
        --out config.json
    """

    stub:
    """
    printf '{"assemblies":[],"tracks":[]}\n' > config.json
    """
}
