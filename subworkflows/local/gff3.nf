// Subworkflow for GFF3 (CAT/Ensembl) → bigBed (bigGenePred)

include { UCSC_GFF3_TO_GENEPRED }          from '../../modules/local/ucsc/gff3ToGenePred.nf'
include { UCSC_GENEPRED_TO_BIGGENEPRED_BED } from '../../modules/local/ucsc/genePredToBigGenePredBed.nf'
include { UCSC_BED_TO_BIGBED_BIGGENEPRED } from '../../modules/local/ucsc/bedToBigBed_bigGenePred.nf'

workflow GFF3_PROCESSING {
    take:
    ch_gff3_files   // channel: [ val(meta), path(gff3) ]
    ch_chrom_sizes  // channel: [ val(genome), path(chrom_sizes) ]

    main:

    // Convert GFF3 -> genePred without touching chrom.sizes
    UCSC_GFF3_TO_GENEPRED(ch_gff3_files)
    ch_gp = UCSC_GFF3_TO_GENEPRED.out.genepred

    UCSC_GENEPRED_TO_BIGGENEPRED_BED(ch_gp)
    ch_bed = UCSC_GENEPRED_TO_BIGGENEPRED_BED.out.biggenepred_bed

    // Join bed with matching chrom.sizes by meta.genome, then add AS file
    ch_sizes_kv2 = ch_chrom_sizes.map { genome, sizes -> [ genome, sizes ] }
    ch_bed_kv = ch_bed.map { meta, bed -> [ meta.genome, [meta, bed] ] }
    // Fan-out chrom.sizes to every BED for the same genome.
    // Do a cross-combine then filter on matching genome keys
    // to avoid join() consuming the single sizes tuple.
    // Result: [meta, bed, sizes] for each (bed, sizes) with genome equality.
    ch_join_sizes = ch_bed_kv
        .combine(ch_sizes_kv2)
        .filter { g_bed, mb, g_sizes, sizes -> g_bed == g_sizes }
        .map { g, mb, g2, sizes -> [ mb[0], mb[1], sizes ] }
    ch_as = Channel.fromPath("${projectDir}/assets/bigGenePred.as")
    // Pair each tuple with the AS file; combine (outer product) is fine because ch_as has one item
    ch_join = ch_join_sizes.combine(ch_as).map { left, as_file -> [ left[0], left[1], left[2], as_file ] }

    UCSC_BED_TO_BIGBED_BIGGENEPRED(ch_join)
    ch_bigbed = UCSC_BED_TO_BIGBED_BIGGENEPRED.out.bigbed

    emit:
    bigbed   = ch_bigbed
}
