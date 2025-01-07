process EXTRACT_SEQUENCES {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::cutadapt=5.0 bioconda::seqkit=2.9.0 conda-forge::pigz=2.8"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_gz)

    output:
    tuple val(sample_id), path("${sample_id}_cut.fastq.gz"), emit: sequences

    script:
    """
    extract_seq_adapter_rc=\$(echo '${params.extract_seq_adapter}' | tr 'ACGTacgt.' 'TGCAtgca.' | rev)

    cutadapt \
        -j $task.cpus \
        -g $params.extract_seq_adapter \
        --minimum-length $params.extract_seq_min_length \
        --maximum-length $params.extract_seq_max_length \
        --discard-untrimmed \
        -o fwd.fastq ${fastq_gz}

    cutadapt \
        -j $task.cpus \
        -g \$extract_seq_adapter_rc \
        --minimum-length $params.extract_seq_min_length \
        --maximum-length $params.extract_seq_max_length \
        --discard-untrimmed \
        -o rev.fastq ${fastq_gz}

    seqkit seq \
        --threads $task.cpus \
        --reverse --complement rev.fastq \
        --out-file rev_rc.fastq

    cat fwd.fastq rev_rc.fastq | pigz -p $task.cpus > ${sample_id}_cut.fastq.gz
    rm fwd.fastq rev.fastq rev_rc.fastq
    """
}