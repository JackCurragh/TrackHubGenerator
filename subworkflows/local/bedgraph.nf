//
// Subworkflow for BEDGRAPH file processing
//

include { CONVERT_BEDGRAPH_TO_BIGWIG }  from '../../modules/local/convert/bedGraphToBigWig.nf'

workflow BEDGRAPH_PROCESSING {
    take:
    ch_bedgraph_files  // channel: [ path(bedgraph_file) ]
    ch_chrom_sizes     // channel: [ path(chrom_sizes) ]

    main:
    ch_versions = Channel.empty()

    // Convert BEDGRAPH to BigWig
    CONVERT_BEDGRAPH_TO_BIGWIG(ch_bedgraph_files, ch_chrom_sizes)
    ch_converted_bigwig = CONVERT_BEDGRAPH_TO_BIGWIG.out.bigwig
    ch_versions = ch_versions.mix(CONVERT_BEDGRAPH_TO_BIGWIG.out.versions)

    emit:
    bigwig   = ch_converted_bigwig  // channel: [ path(bigwig_file) ]
    versions = ch_versions          // channel: [ path(versions.yml) ]
}