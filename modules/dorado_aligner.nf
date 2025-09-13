process DORADO_ALIGNER {
    container 'ontresearch/dorado:latest'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(basecalled), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"), emit: alignment

    script:
    """
    dorado aligner ${reference} ${basecalled} > ${sample_id}.aligned.bam
    """
}
