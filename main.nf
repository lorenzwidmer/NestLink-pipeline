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
    tag "Samtolls fastq on ${basecalled}"

    input:
    path basecalled

    output:
    path "${basecalled}.fastq.gz"

    script:
    """
    samtools fastq $basecalled | gzip > ${basecalled}.fastq.gz
    """
}

/* Workflow */
workflow {
    Channel
        .fromPath(params.data)
        .set { basecalled_channel }
    BamToFastq(basecalled_channel)
}