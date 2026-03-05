//
// Subworkflow for BEDGRAPH file processing (keyed by genome)
//

include { CONVERT_BEDGRAPH_TO_BIGWIG }  from '../../modules/local/convert/bedGraphToBigWig.nf'

workflow BEDGRAPH_PROCESSING {
    take:
    ch_bedgraph_files  // channel: [ val(meta), path(bedgraph_file) ]
    ch_chrom_sizes     // channel: [ val(genome), path(chrom_sizes) ]

    main:

    // Join bedGraph with matching chrom.sizes by meta.genome
    ch_bg_kv   = ch_bedgraph_files.map { meta, bg -> [ meta.genome, [meta, bg] ] }
    ch_sizes_kv = ch_chrom_sizes.map { genome, sizes -> [ genome, sizes ] }
    ch_joined   = ch_bg_kv.join(ch_sizes_kv).map { genome, mb, sizes -> [ mb[0], mb[1], sizes ] }

    // Convert BEDGRAPH to BigWig
    CONVERT_BEDGRAPH_TO_BIGWIG(ch_joined)
    ch_converted_bigwig = CONVERT_BEDGRAPH_TO_BIGWIG.out.bigwig

    emit:
    bigwig   = ch_converted_bigwig  // channel: [ val(meta), path(bigwig_file) ]
}
