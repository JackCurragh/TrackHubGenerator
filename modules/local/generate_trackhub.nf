

    // conda "bioconda::package_name=X.Y.Z"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/package_name:X.Y.Z--hash' :
    //     'biocontainers/package_name:X.Y.Z--hash' }"

process GENERATE_TRACKHUB {
    tag "$hub_name"
    label 'process_low'

    input:
    path(sample_sheet)
    path bigwig
    path bigbed
    val(hub_name)
    val(genome)
    path(output_dir)
    val(sample_regex)
    val(annotation_regex)
    val(email)

    output:
    path "${output_dir}/*"    , emit: trackhub
    path "versions.yml"       , emit: versions

    script:
    def args = task.ext.args ?: ''
    def bigwig_paths = bigwig.collect { "--bigwig '$it'" }.join(' ')
    def bigbed_paths = bigbed.collect { "--bigbed '$it'" }.join(' ')
    """
    TrackHubGenerator.py create \\
        --sample-sheet $sample_sheet \\
        $bigwig_paths \\
        $bigbed_paths \\
        --hub-name "$hub_name" \\
        --genome "$genome" \\
        --output-dir "$output_dir" \\
        --sample-regex "$sample_regex" \\
        --annotation-regex "$annotation_regex" \\
        --email "$email" \\
        $args

    # Version reporting
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TrackHubGenerator: "1.0.0"
        python: \$(python --version | sed 's/Python //')
        trackhub: \$(python -c "import trackhub; print(trackhub.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${output_dir}
    touch ${output_dir}/stub_trackhub.txt
    echo '"${task.process}":' > versions.yml
    echo ' trackhub_creator: "stub_version"' >> versions.yml
    """
}