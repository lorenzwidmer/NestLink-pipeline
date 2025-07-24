process BAM_TO_FASTQ {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::samtools=1.22.1 conda-forge::pigz=2.8"
    tag "${basecalled.baseName}"

    input:
    path(basecalled)

    output:
    tuple val(basecalled.baseName), path("${basecalled.baseName}.fastq.gz"), emit: fastq_gz

    script:
    """
    samtools fastq \
        --threads $task.cpus \
        ${basecalled} \
        | pigz -p $task.cpus > ${basecalled.baseName}.fastq.gz
    """
}