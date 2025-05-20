#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Input parameters
params.bed_files = null
params.bedgraph_files = null
params.bigbed_files = null
params.bigwig_files = null
params.outdir = 'results'
params.genome = null
params.ftp_dir = null

// Import modules
include { VALIDATE_BED } from './modules/local/validate/bed12'
include { VALIDATE_BEDGRAPH } from './modules/local/validate/bedgraph'
include { VALIDATE_BIGWIG } from './modules/local/validate/bigwig'
include { CONVERT_BED_TO_BIGBED } from './modules/local/convert/bedToBigBed'
include { CONVERT_BEDGRAPH_TO_BIGWIG } from './modules/local/convert/bedGraphToBigWig'
include { GET_CHROM_SIZES } from './modules/local/get_chrom_sizes'
include { GENERATE_TRACKHUB } from './modules/local/generate_trackhub'
include { MOVE_TO_FTP } from './modules/local/move_to_ftp'
include { PROCESS_LOG_FILES } from './modules/local/process_log_files'

// Validate inputs
def hasValidInputs = params.bed_files || params.bedgraph_files || params.bigbed_files || params.bigwig_files
if (!hasValidInputs) {
    error "No valid input files provided. Please specify at least one of: bed_files, bedgraph_files, bigbed_files, or bigwig_files"
}

// Main workflow
workflow {
    // Input channels
    ch_bed = params.bed_files ? Channel.fromPath(params.bed_files) : Channel.empty()
    ch_bedgraph = params.bedgraph_files ? Channel.fromPath(params.bedgraph_files) : Channel.empty()
    ch_bigbed = params.bigbed_files ? Channel.fromPath(params.bigbed_files) : Channel.empty()
    ch_bigwig = params.bigwig_files ? Channel.fromPath(params.bigwig_files) : Channel.empty()

    // Validate input files
    ch_validated_bed = VALIDATE_BED(ch_bed)
    ch_validated_bedgraph = VALIDATE_BEDGRAPH(ch_bedgraph)
    ch_validated_bigwig = VALIDATE_BIGWIG(ch_bigwig)

    // Get chromosome sizes
    ch_chrom_sizes = GET_CHROM_SIZES(params.genome)

    // Convert bed to bigBed and bedGraph to bigWig
    ch_converted_bigbed = CONVERT_BED_TO_BIGBED(ch_validated_bed, ch_chrom_sizes)
    ch_converted_bigwig = CONVERT_BEDGRAPH_TO_BIGWIG(ch_validated_bedgraph, ch_chrom_sizes)

    // Collect all bigBed and bigWig files
    ch_all_big_files = ch_converted_bigbed
        .mix(ch_converted_bigwig)
        .mix(ch_validated_bigwig)
        .mix(ch_bigbed)
        .collect()

    // Generate track hub
    ch_trackhub = GENERATE_TRACKHUB(ch_all_big_files)

    // Move to FTP if specified
    if (params.ftp_dir) {
        MOVE_TO_FTP(ch_trackhub, params.ftp_dir)
    }

    // Process log files
    ch_log_files = Channel.fromPath("${params.outdir}/**/*.log")
    PROCESS_LOG_FILES(ch_log_files.collect())
}

// Workflow completion handler
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    
    // Additional log processing if needed
    println "Processing log files..."
    // Add your log processing logic here
}