process MEDAKA_CONSENSUS {
    container 'ontresearch/medaka:sha447c70a639b8bcf17dc49b51e74dfcde6474837b'
    cpus 8
    memory '16 GB'
    time '60m'
    clusterOptions '--gpus=1'
    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(sample_id), path(references), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: consensus
    tuple path("medaka_interference.log"), path("medaka_sequence.log"), emit: log

    script:
    """
    medaka inference \
        --batch 200 --threads 2 \
        --model r1041_e82_400bps_sup_v5.0.0  \
        merged.sorted.bam results.contigs.hdf \
        2> medaka_interference.log

    medaka sequence \
        results.contigs.hdf ${references} ${sample_id}_assembly.fasta \
        2> medaka_sequence.log
    """
}