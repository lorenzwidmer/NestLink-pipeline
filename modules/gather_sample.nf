process GATHER_SAMPLE {
    cpus 1
    memory '1 MB'
    time '10s'
    tag "${basecalled.baseName}"

    input:
    path basecalled

    output:
    tuple val(basecalled.baseName), path(basecalled), emit: calls

    script:
    """
    sleep 1
    """
}
