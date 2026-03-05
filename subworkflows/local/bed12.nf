//
// Subworkflow for BED12 file processing (keyed by genome)
//

include { CONVERT_BED_TO_BIGBED }   from '../../modules/local/convert/bedToBigBed.nf'

workflow BED12_PROCESSING {
    take:
    ch_bed_files     // channel: [ val(meta), path(bed) ]
    ch_chrom_sizes   // channel: [ val(genome), path(chrom_sizes) ]

    main:

    // Join BEDs with their matching chrom.sizes by meta.genome
    ch_bed_kv   = ch_bed_files.map { meta, bed -> [ meta.genome, [meta, bed] ] }
    ch_sizes_kv = ch_chrom_sizes.map { genome, sizes -> [ genome, sizes ] }
    ch_joined   = ch_bed_kv.join(ch_sizes_kv).map { genome, mb, sizes -> [ mb[0], mb[1], sizes ] }

    // Convert BED to BigBed
    CONVERT_BED_TO_BIGBED(ch_joined)
    ch_converted_bigbed = CONVERT_BED_TO_BIGBED.out.bigbed 

    emit:
    bigbed   = ch_converted_bigbed // channel: [ val(meta), path(bigbed_file) ]
}
