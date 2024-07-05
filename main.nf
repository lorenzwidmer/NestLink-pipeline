#!/usr/bin/env nextflow

/* Define parameters */
params.data = "$projectDir/data/*.bam"
params.outdir = "results"

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
    filtlong --min_mean_q 97.5 --min_length 7500 --max_length 7900 ${reads} \
        | gzip > ${baseName}_filtered.fastq.gz
    """
}

process Cutadapt {
    cpus 8
    conda "bioconda::cutadapt=4.9"
    tag "Cutadapt on ${fastq_gz}"

    input:
    path fastq_gz

    output:
    path "*.fastq.gz"

    script:
    """
    cutadapt -j $task.cpus -g atgccatagcatttttatcc...agcctgatacagattaaatc --minimum-length 3839 --maximum-length 3939 --discard-untrimmed -o fwd.fastq.gz ${fastq_gz}
    cutadapt -j $task.cpus -g gatttaatctgtatcaggct...ggataaaaatgctatggcat --minimum-length 3839 --maximum-length 3939 --discard-untrimmed -o rev.fastq.gz ${fastq_gz}
    """
}

/* Workflow */
workflow {
    Channel
        .fromPath(params.data)
        .set { basecalled_ch }
    fastq_gz_ch = BamToFastq(basecalled_ch)
    filtered_ch = FilterReads(fastq_gz_ch)
    Cutadapt(filtered_ch)
}