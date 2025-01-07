process EXTRACT_FLYCODES {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::cutadapt=5.0"
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
}