# NestLink-pipeline
NestLink-pipeline is a pipeline for processing [NestLink libraries](https://www.nature.com/articles/s41592-019-0389-8) sequenced by nanopore sequencing. Reads are binned according to their flycodes (UMIs). Accurate consensus sequences are calculated using Medaka. Variants are called with the pipeline, resulting in a flycode assignment table that links protein variants to their respective set of flycodes.

## Requirements
### Local and cluster execution
- Nextflow ([Installation guide](https://www.nextflow.io/docs/latest/install.html))
- Mamba/ Conda ([https://conda-forge.org/](https://conda-forge.org/))
- mini_align ([mini_align.sh](https://raw.githubusercontent.com/nanoporetech/pomoxis/master/scripts/mini_align) placed in `projectDir/bin/`)
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