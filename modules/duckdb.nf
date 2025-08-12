process DUCKDB {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "conda-forge::python-duckdb=1.3.2"
    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(reads_csv), path(clusters_csv), path(mapped_reads_csv), path(mapped_reads_filtered_csv), path(variants_csv)

    output:
    tuple val(sample_id), path("${sample_id}.db"), emit: database

    script:
    """
    cat << EOF > query.sql

    CREATE TABLE reads AS
    FROM '${sample_id}_reads.csv';

    CREATE TABLE clusters AS
    FROM '${sample_id}_clusters.csv';

    CREATE TABLE mapped_reads AS
    FROM '${sample_id}_mapped_reads.csv';

    CREATE TABLE mapped_reads_filtered AS
    FROM '${sample_id}_mapped_reads_filtered.csv';

    CREATE TABLE variants AS
    FROM '${sample_id}_variants.csv';

    CREATE VIEW barcodes AS
    FROM reads
    JOIN mapped_reads USING(read_id);

    -- Variant summary: barcodes per variant, reads per variant, reads per barcode
    CREATE VIEW variant_summary AS
    SELECT
        "position",
        reference_aa || "position" || variant_aa AS variant,
        COUNT(DISTINCT barcode) AS barcodes,
        COUNT(read_id) AS reads,
        COUNT(read_id) / COUNT(DISTINCT barcode) AS reads_per_barcode
    FROM variants
    JOIN mapped_reads USING (cluster_id)
    WHERE
        variant_type = 'change'
    GROUP BY
        "position", variant
    ORDER BY
        "position" ASC, variant ASC;

    EOF

    duckdb < query.sql ${sample_id}.db
    """
}
