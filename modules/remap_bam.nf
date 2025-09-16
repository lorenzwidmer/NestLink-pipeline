process REMAP_BAM {
    conda "bioconda::dnaio=1.2.3 bioconda::pysam=0.23.3 bioconda::samtools=1.22.1"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(alignment), path(barcode_map)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

    script:
    """
    remap_bam.py \
        --barcode_map ${barcode_map} \
        --input_bam ${alignment} \
        | samtools sort -@ ${task.cpus} --write-index -o ${sample_id}.bam##idx##${sample_id}.bam.bai
    """
}
