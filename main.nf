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
        --min_length $params.filter_min_length \
        --max_length $params.filter_max_length \
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
    extract_gene_adapter_rc=\$(echo '${params.extract_gene_adapter}' | tr 'ACGTacgt.' 'TGCAtgca.' | rev)
    echo \$extract_gene_adapter_rc

    cutadapt \
        -j $task.cpus \
        -g $params.extract_gene_adapter \
        --minimum-length $params.extract_gene_min_length \
        --maximum-length $params.extract_gene_max_length \
        --discard-untrimmed \
        -o fwd.fastq ${fastq_gz}

    cutadapt \
        -j $task.cpus \
        -g \$extract_gene_adapter_rc \
        --minimum-length $params.extract_gene_min_length \
        --maximum-length $params.extract_gene_max_length \
        --discard-untrimmed \
        -o rev.fastq ${fastq_gz}

    seqkit seq \
        --threads $task.cpus \
        --reverse --complement rev.fastq \
        --out-file rev_rc.fastq

    cat fwd.fastq rev_rc.fastq | pigz -p $task.cpus > ${sample_id}_cut.fastq.gz
    rm fwd.fastq rev.fastq rev_rc.fastq
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
        -g $params.extract_flycode_adapter \
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

    publishDir params.outdir, mode: 'copy', pattern: '*.csv'

    input:
    tuple val(sample_id), path(flycodes), path(sequences)
    path(reference)

    output:
    tuple val(sample_id), path("clusters/*.fastq.gz"), path("references/*.fasta"), path("references.fasta"), emit:grouped_reads
    tuple path("flycodes.csv"), path("clusters.csv"), path("mapped_flycodes.csv"), path("mapped_flycodes_filtered.csv"), emit:csv

    script:
    """
    group_by_flycodes.py \
        --flycodes ${flycodes} \
        --sequence ${sequences} \
        --reference_seq ${reference} \
        --reference_flycode ${params.reference_flycode}
    """

    stub:
    """
    mkdir clusters references
    touch clusters/ffffffff-ffff-ffff-ffff-ffffffffffff.fastq.gz references/ffffffff-ffff-ffff-ffff-ffffffffffff.fasta references.fasta
    touch flycodes.csv clusters.csv mapped_flycodes.csv mapped_flycodes_filtered.csv
    """
}

process ALIGN_SEQUENCES {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::minimap2=2.28 bioconda::samtools=1.20"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path("clusters/*"), path("references/*"), path("references.fasta")

    output:
    tuple val(sample_id), path("references.fasta"), path("merged.sorted.bam"), path("merged.sorted.bam.bai"), emit: alignment

    script:
    """
    prepare_alignments.sh $task.cpus

    merge_alignments.py

    rm temp/*.bam
    """

    stub:
    """
    touch merged.sorted.bam merged.sorted.bam.bai
    """
}

process MEDAKA_CONSENSUS {
    container 'ontresearch/medaka:latest'
    cpus 8
    memory '16 GB'
    time '60m'
    clusterOptions '--gpus=1'
    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(sample_id), path(references), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: consensus
    tuple path("medaka_interference.log"), path("medaka_sequence.log"), emit: log

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
    """

    stub:
    """
    touch ${sample_id}_assembly.fasta
    touch medaka_interference.log medaka_sequence.log
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
    def orf2_options = params.orf2_name && params.orf2_pattern ? 
        "--orf2_name ${params.orf2_name} --orf2_pattern ${params.orf2_pattern}" : ''
    """
    variant_calling.py \
        --assembly_path ${assembly} \
        --reference_path ${reference} \
        --experiment_name ${params.experiment_name} \
        --output ${sample_id}_fc.fasta \
        --flycode_pattern ${params.flycode_pattern} \
        --orf1_name ${params.orf1_name} \
        --orf1_pattern ${params.orf1_pattern} \
        ${orf2_options}
    """
    stub:
    """
    touch ${sample_id}_fc.fasta
    """
}

/* Workflows */
workflow nestlink {
    take:
    filter_in_ch
    reference_ch

    main:
    FILTER_READS(filter_in_ch)
    EXTRACT_SEQUENCES(FILTER_READS.out.reads)
    EXTRACT_FLYCODES(EXTRACT_SEQUENCES.out.sequences)
    flycodes_sequences_ch = EXTRACT_FLYCODES.out.flycodes.join(EXTRACT_SEQUENCES.out.sequences)
    GROUP_BY_FLYCODES(flycodes_sequences_ch, reference_ch)
    ALIGN_SEQUENCES(GROUP_BY_FLYCODES.out.grouped_reads)
    MEDAKA_CONSENSUS(ALIGN_SEQUENCES.out.alignment)
    FLYCODE_TABLE(MEDAKA_CONSENSUS.out.consensus, reference_ch)
}

workflow {
    log.info """
    ┌───────────────────────────────────┐
    │ N E S T L I N K   P I P E L I N E │
    │ by Fabian Ackle                   │  
    └───────────────────────────────────┘
    """
    .stripIndent()

    sample_id_ch = Channel.value("barcode05")
    basecalled_ch = Channel.fromPath(params.data)
    filter_in_ch = sample_id_ch.combine(basecalled_ch)
    reference_ch = Channel.fromPath(params.reference)

    nestlink(filter_in_ch, reference_ch)
}