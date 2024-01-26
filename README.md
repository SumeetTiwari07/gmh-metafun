# gmh-metafun (beta)

Nextflow pipeline to perform functional profiling (genefamilies and pathways) of metagenomic short reads.
1. Read cleaning ([fastp](https://github.com/OpenGene/fastp))
2. Summary of Species relative abundance ([metaphlan3](https://huttenhower.sph.harvard.edu/metaphlan3/))
3. Summary of Genefamilies and pathways abundance ([humann3](https://huttenhower.sph.harvard.edu/humann/))
4. Visulaziation of top N species, Genefamiles and pathways.

# Rationale

This is a pipeline for functional as well species abundance profiling of Shotgun Metagenomic reads. This facilitates the estimation of Genefamiles and Pathways' abundance within the samples.

# Dependencies

The pipeline is written in [Nextflow](https://www.nextflow.io/) and its dependencies are available in singularity container.

The YAML file with conda enviromment with required tools is provided in `envs/env.yaml`.

# Databases

# Usage

```
nextflow run humann.nf --reads 'data/*_R{1,2}.fastq.gz' \
    --outdir "nf-humann" -profile nbi
```

Notable options:
* `--uniref [path]`: Path to the uniref database
* `--chocophlan [path]` : Path to the chocophlan database
* `--metaphlandb [path]` : Path to metaphlandb database

Profiles:
* `-profile nbi`: use the default location in the NBI cluster and SLURM scheduler.
* `-profile vmqib`: use the default location of the database in the NBI cluster and local scheduler.

Note: The default location of these database on NBI cluster can be found in the `nextflow.config`
