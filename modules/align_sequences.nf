process ALIGN_SEQUENCES {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::minimap2=2.30 bioconda::samtools=1.22.1"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path("grouped/*"), path("references/*"), path("references.fasta")

    output:
    tuple val(sample_id), path("references.fasta"), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: alignment

    script:
    """
    prepare_alignments.sh -i grouped -r references -o alignments -t ${task.cpus}

    merge_alignments.py -i alignments -o merged.bam

    samtools sort -@ ${task.cpus} merged.bam -o ${sample_id}.bam
    samtools index -@ ${task.cpus} ${sample_id}.bam
    """
}
