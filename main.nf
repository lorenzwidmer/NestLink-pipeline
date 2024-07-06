#!/usr/bin/env nextflow

/* Define parameters */
params.data = "$projectDir/data/*.bam"
params.outdir = "$projectDir/results"

/* Print pipeline info */
log.info """
    =================================
    N E S T L I N K   P I P E L I N E
    =================================
    Data : ${params.data}
    Output dir: ${params.outdir}
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
    path "${basecalled.baseName}.fastq.gz"

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
    path reads

    output:
    path "*_filtered.fastq.gz"

    script:
    def baseName = reads.name.tokenize('.')[0]
    """
    filtlong --min_mean_q 96.84 --min_length 7500 --max_length 7900 ${reads} \
        | gzip > ${baseName}_filtered.fastq.gz
    """
}

process Cutadapt {
    cpus 8
    conda "bioconda::cutadapt=4.9 bioconda::seqkit=2.8.2"
    tag "Cutadapt on ${fastq_gz}"

    input:
    path fastq_gz

    output:
    path "*_cut.fastq.gz"

    script:
    def baseName = fastq_gz.name.tokenize('_')[0]
    """
    cutadapt -j $task.cpus -g atgccatagcatttttatcc...agcctgatacagattaaatc --minimum-length 3839 --maximum-length 3939 --discard-untrimmed -o fwd.fastq ${fastq_gz}
    cutadapt -j $task.cpus -g gatttaatctgtatcaggct...ggataaaaatgctatggcat --minimum-length 3839 --maximum-length 3939 --discard-untrimmed -o rev.fastq ${fastq_gz}

    seqkit seq --threads $task.cpus --reverse --complement rev.fastq --out-file rev_rc.fastq

    cat fwd.fastq rev_rc.fastq > ${baseName}_cut.fastq
    gzip ${baseName}_cut.fastq
    """
}

process NestLink {
    cpus 8
    conda "bioconda::cutadapt=4.9 bioconda::minimap2=2.28 bioconda::samtools=1.20 bioconda::vsearch=2.28"

    input:
    path fastq_gz

    script:
    """
    cutadapt -j $task.cpus -g cagggccccTCAAGA...GGCCAAGGGGGTCAC --error-rate 0.2 --minimum-length 30 --maximum-length 50 --discard-untrimmed --fasta ${fastq_gz} > flycodes.fasta
    mkdir flycode_clusters
    vsearch -cluster_size \
        flycodes.fasta \
        -id 0.90 \
        -sizeout \
        -clusterout_id \
        -centroids flycode_centroids.fasta \
        -clusters flycode_clusters/cluster.fasta
    """
}

/* Workflow */
workflow {
    Channel
        .fromPath(params.data)
        .set { basecalled_ch }
    fastq_gz_ch = BamToFastq(basecalled_ch)
    filtered_ch = FilterReads(fastq_gz_ch)
    cut_ch = Cutadapt(filtered_ch)
    NestLink(cut_ch)
}