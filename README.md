# NestLink-pipeline
NestLink-pipeline is a pipeline for processing [NestLink libraries](https://www.nature.com/articles/s41592-019-0389-8) sequenced by nanopore sequencing. Reads are binned according to their flycodes (UMIs). Accurate consensus sequences are calculated using Medaka. Variants are called with the pipeline, resulting in a flycode assignment table that links protein variants to their respective set of flycodes.

> [!WARNING]
> NestLink-pipeline is still in development. Certain library-specific strings are still hard-coded in `main.nf` and have to be edited before running the pipeline.

## Requirements
- Conda ([https://conda-forge.org/](https://conda-forge.org/))
- Nextflow ([Installation guide](https://www.nextflow.io/docs/latest/install.html))
- Medaka (Note: Medaka is not yet integrated and must be run separately)

## Running the pipeline
1. Clone the repository.
2. Place the basecalled sequencing data and the reference sequence into `projectDir/data/`.
3. Run the first workflow "prepare_data" of the pipeline:
`nextflow run main.nf -entry prepare_data`
4. Generate the consenus sequences using medaka with the data from `projectDir/medaka_input/`, and place the Medaka output `assembly.fasta` into the folder `projectDir/medaka_input/`.
5. Run the second workflow "nestlink" of the pipeline:
`nextflow run main.nf -entry nestlink`

## Generating consensus sequences using Medaka
Example with CUDA and Singularity installed on Ubuntu 20.04.
```bash
singularity run --nv \
    --bind /home/ubuntu/calculation/consensus:/data --pwd /data \
    docker://ontresearch/medaka:latest medaka consensus \
    --batch 200 --threads 2 --model r1041_e82_400bps_sup_v5.0.0  \
    merged.sorted.bam results.contigs.hdf

singularity run --nv \
    --bind /home/ubuntu/calculation/consensus:/data --pwd /data \
    docker://ontresearch/medaka:latest medaka stitch \
    results.contigs.hdf reference.fasta assembly.fasta
```
