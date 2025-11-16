process VARIANT_CALLING {
    conda "bioconda::dnaio=1.2.3 conda-forge::biopython=1.86 conda-forge::polars=1.35.1"
    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy', pattern: '*.txt.gz'

    input:
    tuple val(sample_id), path(assembly), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_variants.csv"), emit: variants_db
    tuple val(sample_id), path("${sample_id}_barcodemap.txt.gz"), emit: enrich2_barcodemap

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
