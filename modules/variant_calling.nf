process VARIANT_CALLING {
    cpus 1
    memory '16 GB'
    time '60m'
    conda "bioconda::dnaio=1.2.3 conda-forge::biopython=1.85 conda-forge::polars=1.26.0"
    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(assembly), path(reference)

    output:
    path "${sample_id}_variants.csv", emit: variants_db
    path "${sample_id}_barcodemap.txt.gz", emit: enrich2_barcodemap

    script:
    """
    variant_calling.py \
        --assembly_path ${assembly} \
        --reference_path ${reference} \
        --sample_id ${sample_id} \
        --barcode_pattern ${params.barcode_pattern.join(' ')} \
        --orf_pattern ${params.orf_pattern.join(' ')} \
        ${params.translate_barcode ? '--translate_barcode' : ''}
        
    gzip ${sample_id}_barcodemap.txt
    """
}