process ALIGN_SEQUENCES {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::minimap2=2.30 bioconda::samtools=1.22.1"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path("clusters/*"), path("references/*"), path("references.fasta")

    output:
    tuple val(sample_id), path("references.fasta"), path("merged.sorted.bam"), path("merged.sorted.bam.bai"), emit: alignment

    script:
    """
    prepare_alignments.sh $task.cpus

    merge_alignments.py

    rm temp/*.bam merged.bam
    """
}