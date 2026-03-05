#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import subworkflows
include { BED12_PROCESSING } from './subworkflows/local/bed12.nf'
include { BEDGRAPH_PROCESSING } from './subworkflows/local/bedgraph.nf'
include { GFF3_PROCESSING } from './subworkflows/local/gff3.nf'

// Import remaining modules
include { GET_CHROM_SIZES_UCSC } from './modules/local/get_chrom_sizes_ucsc.nf'
include { GET_CHROM_SIZES_FASTA } from './modules/local/get_chrom_sizes_fasta.nf'
include { GET_CHROM_SIZES_ASSEMBLY_REPORT } from './modules/local/get_chrom_sizes_assembly_report.nf'
include { GENERATE_TRACKHUB } from './modules/local/generate_trackhub.nf'
include { MOVE_TO_FTP } from './modules/local/move_to_ftp.nf'
include { PROCESS_LOG_FILES } from './modules/local/process_logfiles.nf'
include { AGGREGATE_TRACKHUB } from './modules/local/aggregate_trackhubs.nf'

// Main workflow
workflow {
    // Validate inputs

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
            gff3: filetype == 'gff3'
                return [run, filepath]
        }

    // Access individual file type channels
    ch_bed = ch_files.bed
    ch_bedgraph = ch_files.bedgraph
    ch_bigbed = ch_files.bigbed
    ch_bigwig = ch_files.bigwig
    ch_gff3 = ch_files.gff3

    // Get chromosome sizes
    if (params.chrom_sizes) {
        ch_chrom_sizes = Channel.fromPath(params.chrom_sizes).first()
    } else if (params.assembly_report) {
        def ch_report = Channel.value(file(params.assembly_report))
        def ch_gff_hint = ch_samplesheet
            .filter { meta, ft, p -> ft == 'gff3' }
            .map { meta, ft, p -> p }
            .first()
            .ifEmpty('')
        GET_CHROM_SIZES_ASSEMBLY_REPORT(ch_report, ch_gff_hint)
        ch_chrom_sizes = GET_CHROM_SIZES_ASSEMBLY_REPORT.out.chrom_sizes.first()
    } else if (params.genome) {
        GET_CHROM_SIZES_UCSC(params.genome)
        ch_chrom_sizes = GET_CHROM_SIZES_UCSC.out.chrom_sizes.first()
    } else if (params.genome_fasta) {
        GET_CHROM_SIZES_FASTA(params.genome_fasta)
        ch_chrom_sizes = GET_CHROM_SIZES_FASTA.out.chrom_sizes.first()
    } else {
        error "No chrom sizes provided. Specify one of: --chrom_sizes, --genome_fasta, or --genome."
    }


    // Process BED12 files
    BED12_PROCESSING(ch_bed, ch_chrom_sizes)
    ch_bigbed = ch_bigbed.mix(BED12_PROCESSING.out.bigbed)

    // Process BEDGRAPH files
    BEDGRAPH_PROCESSING(ch_bedgraph, ch_chrom_sizes)
    ch_bigwig = ch_bigwig.mix(BEDGRAPH_PROCESSING.out.bigwig)

    // Process GFF3 annotation files (CAT/Ensembl)
    if (ch_gff3) {
        GFF3_PROCESSING(ch_gff3, ch_chrom_sizes)
        ch_bigbed = ch_bigbed.mix(GFF3_PROCESSING.out.bigbed)
    }

    // Generate track hub - wait for all processing to complete
    ch_trackhub = GENERATE_TRACKHUB(
        ch_bigbed.map { meta, filepath -> filepath }.collect().ifEmpty([]).collect(),
        ch_bigwig.map { meta, filepath -> filepath }.collect().ifEmpty([]).collect(),
        params.hub_name,
        params.genome,
        params.sample_regex,
        params.annotation_regex,
        params.email
    )

    // Optionally aggregate multiple hubs using a JSON manifest of entries
    // Manifest format: [ {"genome": "GCA_...", "trackdb": "/abs/path/to/hub/trackDb.txt"}, ... ]
    if (params.aggregate_name && params.aggregate_manifest) {
        def manifest_path = file(params.aggregate_manifest)
        AGGREGATE_TRACKHUB(
            manifest_path,
            params.aggregate_name,
            params.outdir,
            params.email,
            params.aggregate_short_label ?: params.aggregate_name,
            params.aggregate_long_label ?: "${params.aggregate_name} aggregated hub"
        )
    }

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

// Standalone entry to only build an aggregated multi-assembly hub from a manifest
workflow AGGREGATOR {
    if (!params.aggregate_name || !params.aggregate_manifest) {
        error "AGGREGATOR entry requires --aggregate_name and --aggregate_manifest"
    }

    ch_manifest = Channel.fromPath(params.aggregate_manifest).first()

    AGGREGATE_TRACKHUB(
        ch_manifest,
        params.aggregate_name,
        params.outdir,
        params.email,
        params.aggregate_short_label ?: params.aggregate_name,
        params.aggregate_long_label ?: "${params.aggregate_name} aggregated hub"
    )
}
