#!/usr/bin/env nextflow

/* Processes */
process BAM_TO_FASTQ {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::samtools=1.20 conda-forge::pigz=2.8"
    tag "${basecalled.baseName}"

    input:
    path(basecalled)

    output:
    tuple val(basecalled.baseName), path("${basecalled.baseName}.fastq.gz"), emit: fastq_gz

    script:
    """
    samtools fastq \
        --threads $task.cpus \
        ${basecalled} \
        | pigz -p $task.cpus > ${basecalled.baseName}.fastq.gz
    """

    stub:
    """
    touch ${basecalled.baseName}.fastq.gz
    """
}

process FILTER_READS {
    cpus 1
    memory '4 GB'
    time '60m'
    conda "bioconda::filtlong=0.2.1 conda-forge::pigz=2.8"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), emit: reads

    script:
    """
    filtlong \
        --min_mean_q 96.84 \
        --min_length 7500 --max_length 7900 \
        ${reads} \
        | pigz -p $task.cpus > ${sample_id}_filtered.fastq.gz
    """

    stub:
    """
    touch ${sample_id}_filtered.fastq.gz
    """
}

process EXTRACT_SEQUENCES {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::cutadapt=4.9 bioconda::seqkit=2.8.2 conda-forge::pigz=2.8"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_gz)

    output:
    tuple val(sample_id), path("${sample_id}_cut.fastq.gz"), emit: sequences

    script:
    """
    cutadapt \
        -j $task.cpus \
        -g atgccatagcatttttatcc...agcctgatacagattaaatc \
        --minimum-length 3839 --maximum-length 3939 \
        --discard-untrimmed \
        -o fwd.fastq ${fastq_gz}

    cutadapt \
        -j $task.cpus \
        -g gatttaatctgtatcaggct...ggataaaaatgctatggcat \
        --minimum-length 3839 --maximum-length 3939 \
        --discard-untrimmed \
        -o rev.fastq ${fastq_gz}

    seqkit seq \
        --threads $task.cpus \
        --reverse --complement rev.fastq \
        --out-file rev_rc.fastq

    cat fwd.fastq rev_rc.fastq | pigz -p $task.cpus > ${sample_id}_cut.fastq.gz
    rm fwd.fastq rev_rc.fastq
    """

    stub:
    """
    touch ${sample_id}_cut.fastq.gz
    """
}

process EXTRACT_FLYCODES {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::cutadapt=4.9"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_gz)

    output:
    tuple val(sample_id), path("${sample_id}_flycodes.fasta"), emit: flycodes

    script:
    """
    cutadapt \
        -j $task.cpus \
        -g cagggccccTCAAGA...GGCCAAGGGGGTCAC \
        --error-rate 0.2 \
        --minimum-length 30 --maximum-length 50 \
        --discard-untrimmed \
        --fasta ${fastq_gz} > ${sample_id}_flycodes.fasta
    """

    stub:
    """
    touch ${sample_id}_flycodes.fasta
    """
}

process GROUP_BY_FLYCODES {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::bwa=0.7.18 bioconda::samtools=1.21 bioconda::dnaio=1.2.2 conda-forge::polars=1.17.1 conda-forge::pyarrow=18.1.0 conda-forge::python-duckdb=1.1.3"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(flycodes), path(sequences)
    path(reference)

    output:
    tuple val(sample_id), path("clusters/*.fastq.gz"), path("references/*.fasta"), path("references.fasta"), emit:grouped_reads
    tuple val(sample_id), path("flycodes.csv"), path("clusters.csv"), path("mapped_flycodes.csv"), emit:csv

    script:
    """
    group_by_flycodes.py \
        --flycodes ${flycodes} \
        --sequence ${sequences} \
        --reference_seq ${reference}
    """

    stub:
    """
    touch clusters/ffffffff-ffff-ffff-ffff-ffffffffffff.fastq.gz references/ffffffff-ffff-ffff-ffff-ffffffffffff.fasta references/reference.fasta
    touch flycodes.csv clusters.csv mapped_flycodes.csv
    """
}

process ALIGN_SEQUENCES {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::minimap2=2.28 bioconda::samtools=1.20"
    tag "${sample_id}"

    publishDir "${params.outdir}", mode: 'copy', enabled: workflow.profile == 'standard'

    input:
    tuple val(sample_id), path("clusters/*"), path("references/*"), path("references.fasta")

    output:
    tuple val(sample_id), path("references.fasta"), path("merged.sorted.bam"), path("merged.sorted.bam.bai"), emit: alignment

    script:
    """
    prepare_alignments.sh $task.cpus

    merge_alignments.py
    """

    stub:
    """
    touch reference.fasta merged.sorted.bam merged.sorted.bam.bai
    """
}

process MEDAKA_CONSENSUS {
    container 'ontresearch/medaka:latest'
    cpus 8
    memory '16 GB'
    time '60m'
    clusterOptions '--gpus=1'
    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(references), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: consensus

    script:
    """
    medaka inference \
        --batch 200 --threads 2 \
        --model r1041_e82_400bps_sup_v5.0.0  \
        merged.sorted.bam results.contigs.hdf \
        2> medaka_interference.log

    medaka sequence \
        results.contigs.hdf ${references} ${sample_id}_assembly.fasta \
        2> medaka_sequence.log

    nvidia-smi > nvidia-smi.txt
    """

    stub:
    """
    touch ${sample_id}_assembly.fasta
    """
}

process FLYCODE_TABLE {
    cpus 1
    memory '16 GB'
    time '60m'
    conda "bioconda::dnaio=1.2.1 conda-forge::biopython=1.84"
    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)
    path(reference)

    output:
    path "${sample_id}_fc.fasta", emit: flycode_db

    script:
    """
    flycode_assignment.py \
        --poi "TM287/288_FC" \
        --assembly ${assembly} \
        --file_name ${sample_id}_fc.fasta \
        --reference ${reference}
    """
    stub:
    """
    touch ${sample_id}_fc.fasta
    """
}

/* Workflows */
workflow prepareData {
    take:
    basecalled_ch
    reference_ch

    main:
    BAM_TO_FASTQ(basecalled_ch)
    FILTER_READS(BAM_TO_FASTQ.out.fastq_gz)
    EXTRACT_SEQUENCES(FILTER_READS.out.reads)
    EXTRACT_FLYCODES(EXTRACT_SEQUENCES.out.sequences)
    flycodes_sequences_ch = EXTRACT_FLYCODES.out.flycodes.join(EXTRACT_SEQUENCES.out.sequences)
    GROUP_BY_FLYCODES(flycodes_sequences_ch, reference_ch)
    ALIGN_SEQUENCES(GROUP_BY_FLYCODES.out.grouped_reads)

    emit:
    alignment = ALIGN_SEQUENCES.out.alignment
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

    prepareData(basecalled_ch, reference_ch)

    if (workflow.profile == 'cluster') {
        MEDAKA_CONSENSUS(prepareData.out.alignment)
        FLYCODE_TABLE(MEDAKA_CONSENSUS.out.consensus, reference_ch)
    }
}