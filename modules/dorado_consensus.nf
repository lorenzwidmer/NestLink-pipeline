process DORADO_CONSENSUS {
    cpus 8
    memory '16 GB'
    time '2h'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_polished.fastq"), emit: consensus

    script:
    """
    dorado polish ${bam} ${reference} \
        --qualities \
        --ignore-read-groups \
        --batchsize 250 \
        > ${sample_id}_polished.fastq     
    """
}
