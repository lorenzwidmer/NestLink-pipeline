process VARIANT_CALLING {
    cpus 1
    memory '16 GB'
    time '60m'
    conda "bioconda::dnaio=1.2.2 conda-forge::biopython=1.85 conda-forge::polars=1.17.1"
    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)
    path(reference)

    output:
    path "${sample_id}_fc.fasta", emit: flycode_db

    script:
    """
    variant_calling.py \
        --assembly_path ${assembly} \
        --reference_path ${reference} \
        --experiment_name ${params.experiment_name} \
        --output ${sample_id}_fc.fasta \
        --flycode_pattern ${params.flycode_pattern.join(' ')} \
        --orf_name ${params.orf_name} \
        --orf_pattern ${params.orf_pattern.join(' ')}
    """
}