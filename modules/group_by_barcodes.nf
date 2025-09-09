process GROUP_BY_BARCODES {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::bwa=0.7.19 bioconda::samtools=1.22.1 bioconda::dnaio=1.2.3 conda-forge::polars=1.32.3 conda-forge::pyarrow=21.0.0 conda-forge::python-duckdb=1.3.2"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(barcodes), path(reference)

    output:
    tuple val(sample_id), path("references.fasta"), emit: references
    tuple val(sample_id), path("${sample_id}_mapped_reads_filtered.csv"), emit: barcode_map
    tuple val(sample_id), path("${sample_id}_reads.csv"), path("${sample_id}_clusters.csv"), path("${sample_id}_mapped_reads.csv"), path("${sample_id}_mapped_reads_filtered.csv"), emit: csv

    script:
    """
    group_by_barcodes.py \
        --sample_id ${sample_id} \
        --barcodes ${barcodes} \
        --reference_seq ${reference} \
        --barcode_regex "${params.barcode_regex}" \
        --threads ${task.cpus}
    """
}
