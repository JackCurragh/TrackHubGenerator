//
// Subworkflow for BIGBED file processing
//

include { VALIDATE_BIGBED } from '../../modules/local/validate/bigbed.nf'

workflow BIGBED_PROCESSING {
    take:
    ch_bigbed_files  // channel: [ path(bigbed_file) ]

    main:
    ch_versions = Channel.empty()

    // Validate BigBed files (assuming you have this module)
    ch_validated_bigbed = VALIDATE_BIGBED(ch_bigbed_files)
    ch_versions = ch_versions.mix(VALIDATE_BIGBED.out.versions)

    emit:
    bigbed   = ch_validated_bigbed  // channel: [ path(bigbed_file) ]
    versions = ch_versions          // channel: [ path(versions.yml) ]
}