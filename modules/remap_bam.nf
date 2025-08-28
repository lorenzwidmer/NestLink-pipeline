process REMAP_BAM {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::dnaio=1.2.3 bioconda::pysam=0.23.3 bioconda::samtools=1.22.1"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(alignment), path(barcode_map), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

    script:
    """
    remap_bam.py \
        --barcode_map ${barcode_map} \
        --reference ${reference} \
        --input_bam ${alignment} \
        --output_bam ${sample_id}_unsorted.bam

    samtools sort -@ ${task.cpus} ${sample_id}_unsorted.bam -o ${sample_id}.bam
    samtools index -@ ${task.cpus} ${sample_id}.bam
    """
}
