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

    // Pair each BED with the matching chrom.sizes by genome without consuming sizes
    ch_join_sizes = ch_bed
        .combine(ch_chrom_sizes)
        .filter { meta, bed, genome, sizes -> meta.genome == genome }
        .map    { meta, bed, genome, sizes -> [ meta, bed, sizes ] }

    // Add the AutoSql schema file to each tuple
    ch_as = Channel.fromPath("${projectDir}/assets/bigGenePred.as")
    ch_join = ch_join_sizes
        .combine(ch_as)
        .map { meta, bed, sizes, as_path -> [ meta, bed, sizes, as_path ] }

    UCSC_BED_TO_BIGBED_BIGGENEPRED(ch_join)
    ch_bigbed = UCSC_BED_TO_BIGBED_BIGGENEPRED.out.bigbed

    emit:
    bigbed   = ch_bigbed
}
