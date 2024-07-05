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

process RemoveBarcodes {
    tag "Cutadapt on ${fastq_gz}"

    input:
    path fastq_gz

    output:
    path "*_cut.fastq.gz"

    script:
    """
    remove_barcodes.sh $fastq_gz
    """
}

/* Workflow */
workflow {
    Channel
        .fromPath(params.data)
        .set { basecalled_channel }
    fastq_gz = BamToFastq(basecalled_channel)
    RemoveBarcodes(fastq_gz)
}