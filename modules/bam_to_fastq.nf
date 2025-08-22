process BAM_TO_FASTQ {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::samtools=1.22.1 conda-forge::pigz=2.8"
    tag "${basecalled.baseName}"

    input:
    tuple val(sample_id), path(basecalled)

    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz"), emit: fastq_gz

    script:
    """
    samtools fastq \
        --threads ${task.cpus} \
        ${basecalled} \
        | pigz -p ${task.cpus} > ${sample_id}.fastq.gz
    """
}
