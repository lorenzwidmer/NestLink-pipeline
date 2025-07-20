#!/usr/bin/env nextflow

/* Modules */
include { BAM_TO_FASTQ      } from './modules/bam_to_fastq.nf'
include { FILTER_READS      } from './modules/filter_reads.nf'
include { EXTRACT_SEQUENCES } from './modules/extract_sequences.nf'
include { EXTRACT_FLYCODES  } from './modules/extract_flycodes.nf'
include { GROUP_BY_FLYCODES } from './modules/group_by_flycodes.nf'
include { ALIGN_SEQUENCES   } from './modules/align_sequences.nf'
include { MEDAKA_CONSENSUS  } from './modules/medaka_consensus.nf'
include { VARIANT_CALLING   } from './modules/variant_calling.nf'

/* Workflows */
workflow nestlink {
    take:
    basecalled_ch
    reference_ch

    main:
    BAM_TO_FASTQ(basecalled_ch)
    FILTER_READS(BAM_TO_FASTQ.out.fastq_gz)
    EXTRACT_SEQUENCES(FILTER_READS.out.reads)
    EXTRACT_FLYCODES(EXTRACT_SEQUENCES.out.sequences)
    flycodes_sequences_ch = EXTRACT_FLYCODES.out.flycodes.join(EXTRACT_SEQUENCES.out.sequences)
    GROUP_BY_FLYCODES(flycodes_sequences_ch.combine(reference_ch))
    ALIGN_SEQUENCES(GROUP_BY_FLYCODES.out.grouped_reads)
    MEDAKA_CONSENSUS(ALIGN_SEQUENCES.out.alignment)
    VARIANT_CALLING(MEDAKA_CONSENSUS.out.consensus.combine(reference_ch))
}

workflow {
    log.info """
    ┌───────────────────────────────────┐
    │ N E S T L I N K   P I P E L I N E │
    │ by Fabian Ackle                   │  
    └───────────────────────────────────┘
    """
    .stripIndent()

    basecalled_ch = Channel.fromPath(params.data)
    reference_ch = Channel.fromPath(params.reference)

    nestlink(basecalled_ch, reference_ch)
}