#!/usr/bin/env nextflow

/* Processes */
process BAM_TO_FASTQ {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::samtools=1.20"
    tag "${basecalled.baseName}"

    input:
    path(basecalled)

    output:
    tuple val(basecalled.baseName), path("${basecalled.baseName}.fastq.gz"), emit: fastq_gz

    script:
    """
    samtools fastq --threads $task.cpus ${basecalled} | gzip > ${basecalled.baseName}.fastq.gz
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
    conda "bioconda::filtlong=0.2.1"
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
        | gzip > ${sample_id}_filtered.fastq.gz
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
    conda "bioconda::cutadapt=4.9 bioconda::seqkit=2.8.2"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_gz)

    output:
    tuple val(sample_id), path("${sample_id}_cut.fastq.gz"), emit: sequences

    script:
    """
    cutadapt -j $task.cpus \
        -g atgccatagcatttttatcc...agcctgatacagattaaatc \
        --minimum-length 3839 --maximum-length 3939 \
        --discard-untrimmed \
        -o fwd.fastq ${fastq_gz}

    cutadapt -j $task.cpus \
        -g gatttaatctgtatcaggct...ggataaaaatgctatggcat \
        --minimum-length 3839 --maximum-length 3939 \
        --discard-untrimmed \
        -o rev.fastq ${fastq_gz}

    seqkit seq --threads $task.cpus \
        --reverse --complement rev.fastq \
        --out-file rev_rc.fastq

    cat fwd.fastq rev_rc.fastq | gzip > ${sample_id}_cut.fastq.gz
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
    cutadapt -j $task.cpus \
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

process CLUSTER_FLYCODES {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::vsearch=2.28"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(flycodes)

    output:
    tuple val(sample_id), path("${sample_id}_flycode_centroids.fasta"), path("${sample_id}_flycode_clusters"), emit: clusters

    script:
    """
    mkdir ${sample_id}_flycode_clusters
    vsearch -cluster_size \
        $flycodes \
        -id 0.90 \
        -sizeout \
        -clusterout_id \
        -centroids ${sample_id}_flycode_centroids.fasta \
        -clusters ${sample_id}_flycode_clusters/cluster.fasta
    """

    stub:
    """
    mkdir ${sample_id}_flycode_clusters
    touch ${sample_id}_flycode_centroids.fasta
    """
}

process GROUP_SEQUENCES {
    cpus 1
    memory '4 GB'
    time '60m'
    conda "bioconda::dnaio=1.2.1 conda-forge::pandas=2.2.2"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(centroids), path(clusters)
    tuple val(sample_id2), path(sequences)
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_binned"), emit: binned_reads

    script:
    """
    mkdir ${sample_id}_binned
    group_by_flycodes.py \
        --centroids $centroids \
        --clusters $clusters \
        --sequences $sequences \
        --reference $reference \
        --outdir ${sample_id}_binned
    """

    stub:
    """
    mkdir ${sample_id}_binned
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
    tuple val(sample_id), path(grouped_sequences)
    path(reference)

    output:
    tuple val(sample_id), path("reference_all.fasta"), path("merged.sorted.bam"), path("merged.sorted.bam.bai"), emit: alignment

    script:
    """
    prepare_alignments.sh $reference $grouped_sequences bam $task.cpus
    cp $grouped_sequences/reference.fasta reference_all.fasta
    merge_alignments.py --bam_files bam
    """

    stub:
    """
    touch reference_all.fasta merged.sorted.bam merged.sorted.bam.bai
    """
}

process MEDAKA_CONSENSUS {
    container 'ontresearch/medaka:latest'
    cpus 8
    memory '16 GB'
    time '60m'
    clusterOptions '--gpus=1'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(reference_all), path(bam), path(bai)

    output:
    tuple val(sample_id), path("assembly.fasta"), emit: consensus

    script:
    """
    medaka inference \
    --batch 200 --threads 2 --model r1041_e82_400bps_sup_v5.0.0  \
    merged.sorted.bam results.contigs.hdf \
    > medaka_interference.log

    medaka sequence \
    results.contigs.hdf reference_all.fasta assembly.fasta \
    > medaka_sequence.log

    nvidia-smi > nvidia-smi.txt
    """

    stub:
    """
    touch assembly.fasta
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
    CLUSTER_FLYCODES(EXTRACT_FLYCODES.out.flycodes)
    GROUP_SEQUENCES(CLUSTER_FLYCODES.out.clusters, EXTRACT_SEQUENCES.out.sequences, reference_ch)
    ALIGN_SEQUENCES(GROUP_SEQUENCES.out.binned_reads, reference_ch)

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