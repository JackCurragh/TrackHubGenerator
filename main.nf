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

    // Initialize version tracking
    ch_versions = Channel.empty()

    // Input channels
    ch_samplesheet = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            def meta = [
                id: row.Run,
                filetype: row.FileType
            ]
            [meta, row.FileType, file(row.Path)]
        }

    // Split into separate channels by file type using branch
    ch_files = ch_samplesheet
        .branch { run, filetype, filepath ->
            bed: filetype == 'bed'
                return [run, filepath]
            bedgraph: filetype == 'bedgraph'
                return [run, filepath]
            bigbed: filetype == 'bigbed'
                return [run, filepath]
            bigwig: filetype == 'bigwig'
                return [run, filepath]
        }

    // Access individual file type channels
    ch_bed = ch_files.bed
    ch_bedgraph = ch_files.bedgraph
    ch_bigbed = ch_files.bigbed
    ch_bigwig = ch_files.bigwig

    // Get chromosome sizes
    if (params.genome) {
        GET_CHROM_SIZES_UCSC(params.genome)
        ch_chrom_sizes = GET_CHROM_SIZES_UCSC.out.chrom_sizes.first()
        ch_versions = ch_versions.mix(GET_CHROM_SIZES_UCSC.out.versions)
    } else if (params.genome_fasta) {
        GET_CHROM_SIZES_FASTA(params.genome_fasta)
        ch_chrom_sizes = GET_CHROM_SIZES_FASTA.out.chrom_sizes.first()
        ch_versions = ch_versions.mix(GET_CHROM_SIZES_FASTA.out.versions)
    } else {
        error "No genome or genome fasta provided. Please specify either genome or genome_fasta"
    }


    // Process BED12 files
    BED12_PROCESSING(ch_bed, ch_chrom_sizes)
    ch_bigbed = ch_bigbed.mix(BED12_PROCESSING.out.bigbed)
    ch_versions = ch_versions.mix(BED12_PROCESSING.out.versions)

    // Process BEDGRAPH files
    BEDGRAPH_PROCESSING(ch_bedgraph, ch_chrom_sizes)
    ch_bigwig = ch_bigwig.mix(BEDGRAPH_PROCESSING.out.bigwig)
    ch_versions = ch_versions.mix(BEDGRAPH_PROCESSING.out.versions)

    // Generate track hub
    ch_trackhub = GENERATE_TRACKHUB(
        ch_bigbed.map { meta, filepath -> filepath }.collect().ifEmpty([]),
        [],
        params.hub_name,
        params.genome,
        params.outdir,
        params.sample_regex,
        params.annotation_regex,
        params.email
    )
    ch_versions = ch_versions.mix(GENERATE_TRACKHUB.out.versions)

    // // Move to FTP if specified
    // if (params.ftp_dir) {
    //     MOVE_TO_FTP(ch_trackhub, params.ftp_dir)
    // }

    // // Process log files
    // ch_log_files = Channel.fromPath("${params.outdir}/**/*.log")
    // PROCESS_LOG_FILES(ch_log_files.collect())
}

    // Workflow completion handler
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}