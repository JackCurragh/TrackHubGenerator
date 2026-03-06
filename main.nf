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
    // Expect samplesheet columns: Run,FileType,Path,Genome[,AssemblyReport][,Hub]

    ch_rows = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.Run,
                filetype: row.FileType,
                genome: row.Genome,
                hub: row.containsKey('Hub') ? row.Hub : null,
                assembly_report: row.containsKey('AssemblyReport') ? row.AssemblyReport : ''
            ]
            [ meta, row.FileType, file(row.Path) ]
        }

    // genome -> assembly_report (if provided)
    ch_genome_reports = ch_rows
        .map { meta, ft, p -> [ meta.genome, meta.assembly_report ] }
        .unique()

    // Split by type (meta is carried and includes genome)
    ch_files = ch_rows.branch { meta, filetype, filepath ->
        bed:      filetype == 'bed';       return [meta, filepath]
        bedgraph: filetype == 'bedgraph';  return [meta, filepath]
        bigbed:   filetype == 'bigbed';    return [meta, filepath]
        bigwig:   filetype == 'bigwig';    return [meta, filepath]
        gff3:     filetype == 'gff3';      return [meta, filepath]
    }

    ch_bed      = ch_files.bed
    ch_bedgraph = ch_files.bedgraph
    ch_bigbed   = ch_files.bigbed
    ch_bigwig   = ch_files.bigwig
    ch_gff3     = ch_files.gff3

    // Chrom.sizes per genome via assembly report when available
    ch_from_reports = ch_genome_reports
        .filter { genome, rep -> rep }
        .map    { genome, rep -> [ genome, file(rep) ] }

    // Prepare a GFF hint per genome, defaulting to '' if no GFF is present
    ch_all_genomes = ch_genome_reports.map { g, rep -> g }.unique()
    ch_gff_hint_present = ch_gff3
        .map { meta, p -> [ meta.genome, p ] }
    ch_gff_hint_empty = ch_all_genomes.map { g -> [ g, '' ] }
    ch_gff_hint_all = ch_gff_hint_empty
        .mix(ch_gff_hint_present)
        .groupTuple()
        .map { genome, vals ->
            def v = (vals instanceof List) ? vals.find { it } : vals
            [ genome, v ?: '' ]
        }

    ch_cs_from_report = ch_from_reports.join(ch_gff_hint_all).map { genome, rep, gff -> [ genome, rep, gff ] }
    GET_CHROM_SIZES_ASSEMBLY_REPORT(ch_cs_from_report)
    ch_sizes_report = GET_CHROM_SIZES_ASSEMBLY_REPORT.out.chrom_sizes

    // Fallback: UCSC fetch ONLY for UCSC-style DB names (skip GCA_/GCF_ accessions)
    ch_missing_ucsc = ch_genome_reports
        .filter { g, rep -> !rep && !(g ==~ /^(GCA|GCF)_/ ) }
        .map { g, rep -> g }
    GET_CHROM_SIZES_UCSC(ch_missing_ucsc)
    ch_sizes_ucsc = GET_CHROM_SIZES_UCSC.out.chrom_sizes

    // Unified chrom.sizes stream keyed by genome
    ch_chrom_sizes = ch_sizes_report.mix(ch_sizes_ucsc).unique()

    // Process types with keyed chrom.sizes
    if (ch_bed) {
        BED12_PROCESSING(ch_bed, ch_chrom_sizes)
        ch_bigbed = ch_bigbed.mix(BED12_PROCESSING.out.bigbed)
    }

    if (ch_bedgraph) {
        BEDGRAPH_PROCESSING(ch_bedgraph, ch_chrom_sizes)
        ch_bigwig = ch_bigwig.mix(BEDGRAPH_PROCESSING.out.bigwig)
    }

    if (ch_gff3) {
        GFF3_PROCESSING(ch_gff3, ch_chrom_sizes)
        ch_bigbed = ch_bigbed.mix(GFF3_PROCESSING.out.bigbed)
    }

    // Group tracks into hubs: key = [hubName, genome]
    def hubKey = { meta ->
        def base = meta.hub ?: (params.hub_name ? "${params.hub_name}_${meta.genome}" : meta.genome)
        [ base, meta.genome ]
    }

    ch_bb_kv = ch_bigbed.map { meta, p -> [ hubKey(meta), p ] }
    ch_bw_kv = ch_bigwig.map { meta, p -> [ hubKey(meta), p ] }

    // Combine into one keyed stream, tag each entry by type, then group by key
    ch_comb = ch_bb_kv.map { k, p -> [ k, ['bb', p] ] }.mix( ch_bw_kv.map { k, p -> [ k, ['bw', p] ] } )
    ch_grouped = ch_comb.groupTuple()

    ch_hubs = ch_grouped.map { key, items ->
        def hub_name = key[0]
        def genome = key[1]
        def bbs = items.findAll { it[0] == 'bb' }.collect { it[1] }
        def bws = items.findAll { it[0] == 'bw' }.collect { it[1] }
        [ bbs ?: [], bws ?: [], hub_name, genome, params.sample_regex, params.annotation_regex, params.email ]
    }

    GENERATE_TRACKHUB(ch_hubs)

    // Optionally aggregate multiple hubs using a JSON manifest of entries
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
