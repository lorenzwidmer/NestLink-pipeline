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
    //.baseName
    tag "Samtolls fastq on ${basecalled}"

    input:
    path basecalled

    output:
    path "*.fastq.gz"

    script:
    def fastq_file = basecalled.toString().replace('.bam', '.fastq.gz')
    """
    samtools fastq $basecalled | gzip > $fastq_file
    """
}

process RemoveBarcodes {
    tag "Cutadapt on ${fastq_gz}"
   //publishDir "${params.outdir}", mode: 'copy'

    input:
    path fastq_gz

    output:
    path "*_cut.fastq.gz"

    script:
    def file_name = fastq_gz.toString().replace('.fastq.gz', '')
    """
    remove_barcodes.sh $file_name
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