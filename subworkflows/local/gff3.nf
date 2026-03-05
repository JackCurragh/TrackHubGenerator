// Subworkflow for GFF3 (CAT/Ensembl) → bigBed (bigGenePred)

include { UCSC_GFF3_TO_GENEPRED }          from '../../modules/local/ucsc/gff3ToGenePred.nf'
include { UCSC_GENEPRED_TO_BIGGENEPRED_BED } from '../../modules/local/ucsc/genePredToBigGenePredBed.nf'
include { UCSC_BED_TO_BIGBED_BIGGENEPRED } from '../../modules/local/ucsc/bedToBigBed_bigGenePred.nf'

workflow GFF3_PROCESSING {
    take:
    ch_gff3_files   // channel: [ val(meta), path(gff3) ]
    ch_chrom_sizes  // channel: [ val(genome), path(chrom_sizes) ]

    main:

    // Use chrom sizes twice below; duplicate stream to avoid double-consumption
    ch_chrom_sizes.into { ch_sizes_for_gp; ch_sizes_for_bed }

    // Join with matching chrom.sizes to carry genome context forward (sizes not used here)
    ch_gff_kv   = ch_gff3_files.map { meta, gff -> [ meta.genome, [meta, gff] ] }
    ch_sizes_kv = ch_sizes_for_gp.map { genome, sizes -> [ genome, sizes ] }
    ch_joined   = ch_gff_kv.join(ch_sizes_kv).map { genome, mg, sizes -> [ mg[0], mg[1] ] }

    UCSC_GFF3_TO_GENEPRED(ch_joined)
    ch_gp = UCSC_GFF3_TO_GENEPRED.out.genepred

    UCSC_GENEPRED_TO_BIGGENEPRED_BED(ch_gp)
    ch_bed = UCSC_GENEPRED_TO_BIGGENEPRED_BED.out.biggenepred_bed

    // Join bed with matching chrom.sizes by meta.genome, then add AS file
    ch_sizes_kv2 = ch_sizes_for_bed.map { genome, sizes -> [ genome, sizes ] }
    ch_bed_kv = ch_bed.map { meta, bed -> [ meta.genome, [meta, bed] ] }
    ch_join_sizes = ch_bed_kv.join(ch_sizes_kv2).map { genome, mb, sizes -> [ mb[0], mb[1], sizes ] }
    ch_as = Channel.fromPath("${projectDir}/assets/bigGenePred.as")
    ch_join = ch_join_sizes.cross(ch_as).map { t -> [ t[0], t[1], t[2], t[3] ] }

    UCSC_BED_TO_BIGBED_BIGGENEPRED(ch_join)
    ch_bigbed = UCSC_BED_TO_BIGBED_BIGGENEPRED.out.bigbed

    emit:
    bigbed   = ch_bigbed
}
