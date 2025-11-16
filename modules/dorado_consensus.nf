process DORADO_CONSENSUS {
    container 'ontresearch/dorado:latest'
    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy', pattern: '*.fastq.gz'

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_polished.fastq.gz"), emit: consensus

    script:
    """
    dorado polish ${bam} ${reference} \
        --qualities \
        --ignore-read-groups \
        > ${sample_id}_polished.fastq

    gzip ${sample_id}_polished.fastq   
    """
}
