//
// Subworkflow for BIGWIG file processing
//

include { VALIDATE_BIGWIG } from '../../modules/local/validate/bigwig.nf'

workflow BIGWIG_PROCESSING {
    take:
    ch_bigwig_files  // channel: [ path(bigwig_file) ]

    main:
    ch_versions = Channel.empty()

    // Validate BigWig files
    ch_validated_bigwig = VALIDATE_BIGWIG(ch_bigwig_files)
    ch_versions         = ch_versions.mix(VALIDATE_BIGWIG.out.versions)

    emit:
    bigwig   = ch_validated_bigwig  // channel: [ path(bigwig_file) ]
    versions = ch_versions          // channel: [ path(versions.yml) ]
}