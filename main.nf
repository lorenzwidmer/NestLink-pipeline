#!/usr/bin/env nextflow

/* Modules */
include { GATHER_SAMPLE     } from './modules/gather_sample.nf'
include { BAM_TO_FASTQ      } from './modules/bam_to_fastq.nf'
include { DORADO_ALIGNER    } from './modules/dorado_aligner.nf'
include { FILTER_READS      } from './modules/filter_reads.nf'
include { EXTRACT_BARCODES  } from './modules/extract_barcodes.nf'
include { GROUP_BY_BARCODES } from './modules/group_by_barcodes.nf'
include { ALIGN_SEQUENCES   } from './modules/align_sequences.nf'
include { MEDAKA_CONSENSUS  } from './modules/medaka_consensus.nf'
include { VARIANT_CALLING   } from './modules/variant_calling.nf'
include { DUCKDB            } from './modules/duckdb.nf'

/* Workflows */
workflow {
    log.info(
        """
        ┌─────────────────────────────────────────────┐
        │ O N T   C O N S E N S U S   P I P E L I N E │
        │ by Fabian Ackle                             │  
        └─────────────────────────────────────────────┘
        """.stripIndent()
    )

    basecalled_ch = Channel.fromPath(params.data)
    reference_ch = Channel.fromPath(params.reference)

    consensus(basecalled_ch, reference_ch)
}

workflow consensus {
    take:
    basecalled_ch
    reference_ch

    main:
    GATHER_SAMPLE(basecalled_ch)
    BAM_TO_FASTQ(GATHER_SAMPLE.out.calls)
    DORADO_ALIGNER(GATHER_SAMPLE.out.calls.combine(reference_ch))
    FILTER_READS(BAM_TO_FASTQ.out.fastq_gz)
    EXTRACT_BARCODES(FILTER_READS.out.reads)
    barcodes_sequences_ch = EXTRACT_BARCODES.out.barcodes.join(FILTER_READS.out.reads)
    GROUP_BY_BARCODES(barcodes_sequences_ch.combine(reference_ch))
    ALIGN_SEQUENCES(GROUP_BY_BARCODES.out.grouped_reads)
    MEDAKA_CONSENSUS(ALIGN_SEQUENCES.out.alignment)
    VARIANT_CALLING(MEDAKA_CONSENSUS.out.consensus.combine(reference_ch))
    csv_ch = GROUP_BY_BARCODES.out.csv.join(VARIANT_CALLING.out.variants_db)
    DUCKDB(csv_ch)
}
