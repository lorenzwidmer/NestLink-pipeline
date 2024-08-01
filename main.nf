#!/usr/bin/env nextflow

/* Define parameters */
params.data = "$projectDir/data/*.bam"
params.reference = "$projectDir/data/reference.fasta"
params.medeka_in = "$projectDir/medaka_input"
params.medeka_out = "$projectDir/medaka_output/*.fasta"
params.outdir = "$projectDir/results"

/* Print pipeline info */
log.info """
    =================================
    N E S T L I N K   P I P E L I N E
    =================================
    """
    .stripIndent()

/* Processes */
process BamToFastq {
    cpus 8
    conda "bioconda::samtools=1.20"
    tag "Samtools fastq on ${basecalled}"

    input:
    path basecalled

    output:
    tuple val(basecalled.baseName), path("${basecalled.baseName}.fastq.gz")

    script:
    """
    samtools fastq --threads $task.cpus $basecalled | gzip > "${basecalled.baseName}.fastq.gz"
    """
}

process FilterReads {
    cpus 1
    conda "bioconda::filtlong=0.2.1"
    tag "Filtlong on ${reads}"

    input:
    tuple val(sampleName), path(reads)

    output:
    tuple val(sampleName), path("${sampleName}_filtered.fastq.gz")

    script:
    """
    filtlong --min_mean_q 96.84 \
        --min_length 7500 --max_length 7900 \
        ${reads} | gzip > ${sampleName}_filtered.fastq.gz
    """
}

process ExtractSequences {
    cpus 8
    conda "bioconda::cutadapt=4.9 bioconda::seqkit=2.8.2"
    tag "Cutadapt on ${fastq_gz}"

    input:
    tuple val(sampleName), path(fastq_gz)

    output:
    tuple val(sampleName), path("${sampleName}_cut.fastq.gz")

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

    cat fwd.fastq rev_rc.fastq > ${sampleName}_cut.fastq
    gzip ${sampleName}_cut.fastq
    rm *.fastq
    """
}

process ExtractFlycodes {
    cpus 8
    conda "bioconda::cutadapt=4.9"
    tag "Cutadapt on $fastq_gz"

    input:
    tuple val(sampleName), path(fastq_gz)

    output:
    tuple val(sampleName), path("${sampleName}_flycodes.fasta")

    script:
    """
    cutadapt -j $task.cpus \
        -g cagggccccTCAAGA...GGCCAAGGGGGTCAC \
        --error-rate 0.2 \
        --minimum-length 30 --maximum-length 50 \
        --discard-untrimmed \
        --fasta ${fastq_gz} > ${sampleName}_flycodes.fasta
    """
}

process ClusterFlycodes {
    cpus 8
    conda "bioconda::vsearch=2.28"
    tag "Vsearch on $flycodes"

    input:
    tuple val(sampleName), path(flycodes)

    output:
    tuple path("${sampleName}_flycode_centroids.fasta"), path("${sampleName}_flycode_clusters")

    script:
    """
    mkdir ${sampleName}_flycode_clusters
    vsearch -cluster_size \
        $flycodes \
        -id 0.90 \
        -sizeout \
        -clusterout_id \
        -centroids ${sampleName}_flycode_centroids.fasta \
        -clusters ${sampleName}_flycode_clusters/cluster.fasta
    """
}

process GroupSequences {
    cpus 1
    conda "bioconda::dnaio=1.2.1"
    tag "group_by_flycodes.py on $sequences with $centroids"

    input:
    tuple path(centroids), path(clusters)
    tuple val(sampleName), path(sequences)

    output:
    tuple val(sampleName), path("${sampleName}_binned")

    script:
    """
    mkdir ${sampleName}_binned
    group_by_flycodes.py \
        --centroids $centroids \
        --clusters $clusters \
        --sequences $sequences \
        --outdir ${sampleName}_binned
    """
}

process AlignSequences {
    cpus 8
    conda "bioconda::minimap2=2.28 bioconda::samtools=1.20"
    tag "mini_align on $grouped_sequences"

    publishDir "${params.medeka_in}/${sampleName}", mode: 'copy'

    input:
    tuple val(sampleName), path(grouped_sequences), path(reference)

    output:
    path "reference_all.fasta"
    path "merged.sorted.bam"
    path "merged.sorted.bam.bai"

    script:
    """
    prepare_alignments.sh $reference $grouped_sequences bam $task.cpus
    merge_alignments.py --bam_files bam
    """
}

process makeFlycodeTable {
    cpus 1
    conda "bioconda::dnaio=1.2.1 conda-forge::biopython=1.84"
    tag "flycode_assignment.py on $assembly"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple path(assembly), path(reference)

    output:
    path "${assembly.baseName}_fc.fasta"

    script:
    """
    flycode_assignment.py \
        --poi "TM287/288_FC" \
        --assembly $assembly \
        --file_name ${assembly.baseName}_fc.fasta \
        --reference $reference
    """
}

/* Workflow */
workflow prepare_data {
    Channel
        .fromPath(params.data)
        .set { basecalled_ch }
    Channel
        .fromPath(params.reference)
        .set { reference_ch }
    fastq_gz_ch = BamToFastq(basecalled_ch)
    filtered_ch = FilterReads(fastq_gz_ch)
    sequences_ch = ExtractSequences(filtered_ch)
    flycodes_ch = ExtractFlycodes(sequences_ch)
    clusters_ch = ClusterFlycodes(flycodes_ch)
    group_ch = GroupSequences(clusters_ch, sequences_ch)
    align_inp_ch = group_ch.combine(reference_ch)
    /*align_inp_ch.view()*/
    AlignSequences(align_inp_ch)
}

workflow nestlink {
    Channel
        .fromPath(params.medeka_out)
        .set { consensus_ch }
    Channel
        .fromPath(params.reference)
        .set { reference_ch }
    nestlink_inp_ch = consensus_ch.combine(reference_ch)
    makeFlycodeTable(nestlink_inp_ch)
}