

process PROCESS_NAME {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::package_name=X.Y.Z"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/package_name:X.Y.Z--hash' :
        'biocontainers/package_name:X.Y.Z--hash' }"

    input:
    path trackhub_dir
    path ftp_destination

    output:
    path "*.log" , emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    become genebuild 
    
    rsync -avz ${trackhub_dir} ${ftp_destination}
    """
}