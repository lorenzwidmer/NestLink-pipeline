process FLYCODE_TABLE {
    cpus 1
    memory '16 GB'
    time '60m'
    conda "bioconda::dnaio=1.2.2 conda-forge::biopython=1.84"
    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)
    path(reference)

    output:
    path "${sample_id}_fc.fasta", emit: flycode_db

    script:
    def orf2_options = params.orf2_name && params.orf2_pattern ? 
        "--orf2_name ${params.orf2_name} --orf2_pattern ${params.orf2_pattern}" : ''
    """
    variant_calling.py \
        --assembly_path ${assembly} \
        --reference_path ${reference} \
        --experiment_name ${params.experiment_name} \
        --output ${sample_id}_fc.fasta \
        --flycode_pattern ${params.flycode_pattern} \
        --orf1_name ${params.orf1_name} \
        --orf1_pattern ${params.orf1_pattern} \
        ${orf2_options}
    """
}