// Subworkflow for GFF3 (CAT/Ensembl) → bigBed (bigGenePred)

include { UCSC_GFF3_TO_GENEPRED }          from '../../modules/local/ucsc/gff3ToGenePred.nf'
include { UCSC_GENEPRED_TO_BIGGENEPRED_BED } from '../../modules/local/ucsc/genePredToBigGenePredBed.nf'
include { UCSC_BED_TO_BIGBED_BIGGENEPRED } from '../../modules/local/ucsc/bedToBigBed_bigGenePred.nf'

workflow GFF3_PROCESSING {
    take:
    ch_gff3_files   // channel: [ val(meta), path(gff3) ] or [ path(gff3) ]
    ch_chrom_sizes  // channel: [ path(chrom_sizes) ]

    main:
    ch_versions = Channel.empty()

    // Normalize inputs to (meta, path)
    ch_norm = ch_gff3_files.map { item ->
        if (item instanceof Path) {
            def meta = [ id: item.baseName ]
            return [ meta, item ]
        } else if (item instanceof Tuple && item.size()==2) {
            return item
        } else {
            def (meta, p) = item
            return [ meta, p ]
        }
    }

    UCSC_GFF3_TO_GENEPRED(ch_norm)
    ch_gp = UCSC_GFF3_TO_GENEPRED.out.genepred
    ch_versions = ch_versions.mix(UCSC_GFF3_TO_GENEPRED.out.versions)

    UCSC_GENEPRED_TO_BIGGENEPRED_BED(ch_gp)
    ch_bed = UCSC_GENEPRED_TO_BIGGENEPRED_BED.out.biggenepred_bed
    ch_versions = ch_versions.mix(UCSC_GENEPRED_TO_BIGGENEPRED_BED.out.versions)

    // Broadcast constants: chrom.sizes and AS file (Cartesian join)
    ch_sizes = ch_chrom_sizes
    ch_as = Channel.fromPath("${projectDir}/assets/bigGenePred.as")

    ch_join = ch_bed.cross(ch_sizes).cross(ch_as).map { tuple1 ->
        def t1 = tuple1[0]      // [ [meta, bed], sizes ]
        def meta_bed = t1[0]    // [meta, bed]
        def sizes = t1[1]
        def asf = tuple1[1]
        [ meta_bed[0], meta_bed[1], sizes, asf ]
    }

    UCSC_BED_TO_BIGBED_BIGGENEPRED(ch_join)
    ch_bigbed = UCSC_BED_TO_BIGBED_BIGGENEPRED.out.bigbed
    ch_versions = ch_versions.mix(UCSC_BED_TO_BIGBED_BIGGENEPRED.out.versions)

    emit:
    bigbed   = ch_bigbed
    versions = ch_versions
}
