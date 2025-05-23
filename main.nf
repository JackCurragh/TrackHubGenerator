#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import subworkflows
include { BED12_PROCESSING } from './subworkflows/local/bed12.nf'
include { BEDGRAPH_PROCESSING } from './subworkflows/local/bedgraph.nf'

// Import remaining modules
include { GET_CHROM_SIZES_UCSC } from './modules/local/get_chrom_sizes_ucsc.nf'
include { GET_CHROM_SIZES_FASTA } from './modules/local/get_chrom_sizes_fasta.nf'
include { GENERATE_TRACKHUB } from './modules/local/generate_trackhub.nf'
include { MOVE_TO_FTP } from './modules/local/move_to_ftp.nf'
include { PROCESS_LOG_FILES } from './modules/local/process_logfiles.nf'

// Main workflow
workflow {
    // Validate inputs
    def hasValidInputs = params.bed_files || params.bedgraph_files || params.bigbed_files || params.bigwig_files
    if (!hasValidInputs) {
        error "No valid input files provided. Please specify at least one of: bed_files, bedgraph_files, bigbed_files, or bigwig_files"
    }

    // Initialize version tracking
    ch_versions = Channel.empty()

    // Input channels
    ch_bed = params.bed_files ? Channel.fromPath(params.bed_files) : Channel.empty()
    ch_bedgraph = params.bedgraph_files ? Channel.fromPath(params.bedgraph_files) : Channel.empty()
    ch_bigbed = params.bigbed_files ? Channel.fromPath(params.bigbed_files) : Channel.empty()
    ch_bigwig = params.bigwig_files ? Channel.fromPath(params.bigwig_files) : Channel.empty()

    // Get chromosome sizes
    if (params.genome) {
        ch_chrom_sizes = GET_CHROM_SIZES_UCSC(params.genome)
        ch_versions = ch_versions.mix(GET_CHROM_SIZES_UCSC.out.versions)
    } else if (params.genome_fasta) {
        ch_chrom_sizes = GET_CHROM_SIZES_FASTA(params.genome_fasta)
        ch_versions = ch_versions.mix(GET_CHROM_SIZES_FASTA.out.versions)
    } else {
        error "No genome or genome fasta provided. Please specify either genome or genome_fasta"
    }

    // Process BED12 files
    // BED12_PROCESSING(ch_bed, ch_chrom_sizes)
    // ch_bigbed = ch_bigbed.mix(BED12_PROCESSING.out.bigbed)
    // ch_versions = ch_versions.mix(BED12_PROCESSING.out.versions)

    // // Process BEDGRAPH files
    // BEDGRAPH_PROCESSING(ch_bedgraph, ch_chrom_sizes)
    // ch_bigwig = ch_bigwig.mix(BEDGRAPH_PROCESSING.out.bigwig)
    // ch_versions = ch_versions.mix(BEDGRAPH_PROCESSING.out.versions)

    // // Generate track hub
    // ch_trackhub = GENERATE_TRACKHUB(
    //     params.sample_sheet,
    //     ch_bigbed.collect(),
    //     ch_bigwig.collect(),
    //     params.hub_name,
    //     params.genome,
    //     params.outdir,
    //     params.sample_regex,
    //     params.annotation_regex,
    //     params.email
    // )
    // ch_versions = ch_versions.mix(GENERATE_TRACKHUB.out.versions)

    // // Move to FTP if specified
    // if (params.ftp_dir) {
    //     MOVE_TO_FTP(ch_trackhub, params.ftp_dir)
    // }

    // // Process log files
    // ch_log_files = Channel.fromPath("${params.outdir}/**/*.log")
    // PROCESS_LOG_FILES(ch_log_files.collect())

    // Workflow completion handler
    workflow.onComplete {
        log.info "Pipeline completed at: $workflow.complete"
        log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    }
}