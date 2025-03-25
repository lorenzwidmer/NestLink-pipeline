# NestLink-pipeline
NestLink-pipeline is a pipeline for processing [NestLink libraries](https://www.nature.com/articles/s41592-019-0389-8) sequenced by nanopore sequencing. Reads are binned according to their flycodes (UMIs). Accurate consensus sequences are calculated using Medaka. Variants are called with the pipeline, resulting in a flycode assignment table that links protein variants to their respective set of flycodes.

## Requirements
### Local and cluster execution
- Nextflow ([Installation guide](https://www.nextflow.io/docs/latest/install.html)), on the cluster it has to be installed in a mamba/ conda environment called `nextflow`. 
- Mamba/ Conda ([https://conda-forge.org/](https://conda-forge.org/))
- mini_align ([mini_align.sh](https://raw.githubusercontent.com/nanoporetech/pomoxis/master/scripts/mini_align) placed in `./bin/`)
### Local execution only
- Podman ([https://podman.io/](https://podman.io/))
### Cluster execution only
- Slurm workflow manager
- Singularity

## Running the pipeline on the s3it cluster
1. Clone the repository.
2. Edit the params.json file, specify the nanopore reads (bam) and reference sequence.
3. Run the pipeline:
`sbatch run_NL-pipeline.slurm`

## Running the pipeline locally
1. Prepare the pipeline as described above.
2. Run the pipeline:
`bash run_NL-pipeline.sh`

## Parameters
| Parameter                 | Type                 | Description                                 |
|---------------------------|----------------------|---------------------------------------------|
| `data`                    | String               | Path to input BAM file.                     |
| `reference`               | String               | Path to reference FASTA file.               |
| `filter_min_length`       | Integer              | Read filtering minimum length threshold.    |
| `filter_max_length`       | Integer              | Read filtering maximum length threshold.    |
| `extract_seq_adapter`     | String               | Linked adapter for sequence trimming.       |
| `extract_seq_min_length`  | Integer              | Sequence trimming minimum length threshold. |
| `extract_seq_max_length`  | Integer              | Sequence trimming minimum length threshold. |
| `extract_flycode_adapter` | String               | Linked adapter for flycode extraction.      |
| `medaka_dorado_model`     | String               | Dorado model used for basecalling.          |
| `flycode_pattern`         | List(String, String) | Sequences flanking flyodes.                 |
| `orf_pattern`             | List(String, String) | Sequences flanking ORF.                     |
| `outdir`                  | String               | Output directory for results.               |
