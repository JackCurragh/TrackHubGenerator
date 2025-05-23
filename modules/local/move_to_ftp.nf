

process MOVE_TO_FTP {
    label 'process_low'

    input:
    path trackhub_dir
    path ftp_destination

    output:
    path "*.log" , emit: log

    script:
    """
    become genebuild 
    
    rsync -avz ${trackhub_dir} ${ftp_destination}
    """
}