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
include { VALIDATE_BED } from './modules/local/validate/bed12.nf'
include { VALIDATE_BEDGRAPH } from './modules/local/validate/bedgraph.nf'
include { VALIDATE_BIGWIG } from './modules/local/validate/bigwig.nf'
include { CONVERT_BED_TO_BIGBED } from './modules/local/convert/bedToBigBed.nf'
include { CONVERT_BEDGRAPH_TO_BIGWIG } from './modules/local/convert/bedGraphToBigWig'
include { GET_CHROM_SIZES_UCSC } from './modules/local/get_chrom_sizes_ucsc'
include { GET_CHROM_SIZES_FASTA } from './modules/local/get_chrom_sizes_fasta'
include { GENERATE_TRACKHUB } from './modules/local/generate_trackhub'
include { MOVE_TO_FTP } from './modules/local/move_to_ftp'
include { PROCESS_LOG_FILES } from './modules/local/process_logfiles'



// Main workflow
workflow {
    // Validate inputs
    def hasValidInputs = params.bed_files || params.bedgraph_files || params.bigbed_files || params.bigwig_files
    if (!hasValidInputs) {
        error "No valid input files provided. Please specify at least one of: bed_files, bedgraph_files, bigbed_files, or bigwig_files"
    }
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
    if (params.genome) {
            ch_chrom_sizes = GET_CHROM_SIZES_UCSC(params.genome)
    } else if (params.genome_fasta) {
            ch_chrom_sizes = GET_CHROM_SIZES_FASTA(params.genome_fasta)
    } else {
            error "No genome or genome fasta provided. Please specify either genome or genome_fasta"
    }

    // Convert bed to bigBed and bedGraph to bigWig
    ch_converted_bigbed = CONVERT_BED_TO_BIGBED(ch_validated_bed, ch_chrom_sizes)
    ch_converted_bigwig = CONVERT_BEDGRAPH_TO_BIGWIG(ch_validated_bedgraph, ch_chrom_sizes)

    // Generate track hub
    ch_trackhub = GENERATE_TRACKHUB(
        params.sample_sheet,
        ch_converted_bigbed.collect(),
        ch_converted_bigwig.collect(),
        params.hub_name,
        params.genome,
        params.outdir,
        params.sample_regex,
        params.annotation_regex,
        params.email,
        )

    // Move to FTP if specified
    if (params.ftp_dir) {
        MOVE_TO_FTP(ch_trackhub, params.ftp_dir)
    }

    // Process log files
    ch_log_files = Channel.fromPath("${params.outdir}/**/*.log")
    PROCESS_LOG_FILES(ch_log_files.collect())

    // Workflow completion handler
    workflow.onComplete {
        log.info "Pipeline completed at: $workflow.complete"
        log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    }
}

