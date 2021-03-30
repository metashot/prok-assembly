#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { deinterlive; trim_adapters; remove_contaminants; quality_filter;
    raw_reads_stats; clean_reads_stats; statswrapper } from './modules/bbtools'

workflow {
    Channel
        .fromFilePairs( params.reads, size: (params.single_end || params.interleaved) ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}." }
        .set { reads_ch }

    if (params.interleaved) {
        deinterlive(reads_ch)
        deint_ch = deinterlive.out.reads
    } else {
        deint_ch = reads_ch
    }

    raw_reads_stats(deint_ch)

    if (!params.skip_cleaning) {
        trim_adapters(deint_ch)
        remove_contaminants(trim_adapters.out.reads)
        clean_ch = quality_filter(remove_contaminants.out.reads)
        clean_reads_stats(qual_ch)
    } else {
        clean_ch = reads_ch
    }

    if (!params.skip_spades) {
        spades(qual_ch)
        scaffolds_ch = spades.out.scaffolds
        statswrapper(scaffolds_ch.map { row -> row[1] }.collect())
    }

    if (!params.skip_plasmidspades) {
        plasmidspades(qual_ch)

        if (params.viralverify_db == 'none') {
            viralverify_db_download()
            viralverify_db = viralverify_db_download.out.viralverify_db
        }
        else {
            viralverify_db = file(params.viralverify_db, checkIfExists: true)
        }

        viralverify(spades.out.scaffolds)
    }
}
