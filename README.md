# NestLink-pipeline

NestLink-pipeline is a Nextflow pipeline for processing barcoded protein variant libraries and [NestLink libraries](https://www.nature.com/articles/s41592-019-0389-8) sequenced by nanopore sequencing.
Reads are binned according to their barcodes or flycodes (UMIs).
Accurate consensus sequences are calculated using Dorado polish.
Finally, variants are called with the pipeline, linking barcodes or flycodes with their respective protein variants.

## Requirements
### Local and cluster execution
- Nextflow ([Installation guide](https://www.nextflow.io/docs/latest/install.html)), on the cluster it has to be installed in a mamba/ conda environment called `nextflow`. 
- Mamba/ Conda ([https://conda-forge.org/](https://conda-forge.org/))
### Local execution only
- Dorado ([Installation guide](https://software-docs.nanoporetech.com/dorado/latest/#installation)).
### Cluster execution only
- Slurm workflow manager
- Singularity

## Running the pipeline

1. Clone the repository with `git clone https://github.com/fabianackle/NestLink-pipeline.git`.
2. Create a `params.json` file with the parameters listed below, specify the nanopore reads (BAM) and reference sequence, see the examples contained in this repo.
3. Run the pipeline with either `./run_NL-pipeline.sh` for local execution or on a cluster with `sbatch run_NL-pipeline.slurm`.

## Parameters

| Parameter                 | Type                 | Description                                                                         |
|---------------------------|----------------------|-------------------------------------------------------------------------------------|
| `data`                    | String               | Path to input BAM file(s).†                                                         |
| `reference`               | String               | Path to reference FASTA file.                                                       |
| `filter_quality`          | Float                | Minimum mean read quality threshold.                                                |
| `filter_min_length`       | Integer              | Read filtering minimum length threshold.                                            |
| `filter_max_length`       | Integer              | Read filtering maximum length threshold.                                            |
| `extract_barcode_adapter` | String               | Linked adapter for barcode extraction.                                              |
| `barcode_regex`           | String               | Regular expression matching the barcode.                                            |
| `barcode_min_coverage`    | Integer              | The minimal amount a barcode has to be seen to be considered a high-quality barcode.|
| `barcode_pattern`         | List(String, String) | Sequences flanking the barcode.                                                     |
| `orf_pattern`             | List(String, String) | Sequences flanking ORF.                                                             |
| `translate_barcode`       | Boolean              | Translates barcode, used with flycodes.                                             |
| `outdir`                  | String               | Output directory for results.                                                       |

† for multiple BAM files use `*`, e.g. `data/barcode*.bam`.
