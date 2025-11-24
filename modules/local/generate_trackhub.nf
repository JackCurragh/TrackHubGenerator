process GENERATE_TRACKHUB {
    tag "$hub_name"
    label 'process_low'

    container "community.wave.seqera.io/library/pip_trackhub:b1b9686e5cada428"

    publishDir "${params.outdir}/trackhubs", mode: 'copy'

    input:
    path(bigbed)
    path(bigwig)
    val(hub_name)
    val(genome)
    val(sample_regex)
    val(annotation_regex)
    val(email)

    output:
    path "trackhub_output/**"              , emit: trackhub
    path "versions.yml"                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    def bigwig_paths = bigwig ? "--bigwig " + bigwig.collect { "'$it'" }.join(' ') : ''
    def bigbed_paths = bigbed ? "--bigbed " + bigbed.collect { "'$it'" }.join(' ') : ''
    def annotation_regex_param = annotation_regex ? "--annotation-regex '$annotation_regex'" : ''
    def sample_regex_param = sample_regex ? "--sample-regex '$sample_regex'" : ''
    """
    TrackHubGenerator.py create \\
        $bigwig_paths \\
        $bigbed_paths \\
        --hub-name "$hub_name" \\
        --genome "$genome" \\
        --output-dir trackhub_output \\
        $sample_regex_param \\
        $annotation_regex_param \\
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
    mkdir -p trackhub_output/${hub_name}
    touch trackhub_output/${hub_name}/stub_trackhub.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TrackHubGenerator: "stub_version"
    END_VERSIONS
    """
}