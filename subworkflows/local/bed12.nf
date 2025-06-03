//
// Subworkflow for BED12 file processing
//

include { CONVERT_BED_TO_BIGBED }   from '../../modules/local/convert/bedToBigBed.nf'

workflow BED12_PROCESSING {
    take:
    ch_bed_files     // channel: [ path(bed_file) ]
    ch_chrom_sizes   // channel: [ path(chrom_sizes) ]

    main:
    ch_versions = Channel.empty()

    // Convert BED to BigBed
    CONVERT_BED_TO_BIGBED(ch_bed_files, ch_chrom_sizes)
    ch_converted_bigbed = CONVERT_BED_TO_BIGBED.out.bigbed 
    ch_versions = ch_versions.mix(CONVERT_BED_TO_BIGBED.out.versions)

    emit:
    bigbed   = ch_converted_bigbed // channel: [ path(bigbed_file) ]
    versions = ch_versions          // channel: [ path(versions.yml) ]
}