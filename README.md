# Humann pipeline (beta)

Nextflow pipeline to perform functional profiling (genefamilies and pathways) of metagenomic short reads.
1. Read cleaning (fastp)
2. Summary of Species relative abundance (metaphlan3)
3. Summary of Genefamilies and pathways abundance (humann3)
4. Visulaziation of top N species, Genefamiles and pathways.


How to run the pipeline:

1. Install Nextflow (>=21.04.0)
2. Create environment and activate it (`conda env create -n gmh-humann --file envs/env.yaml`)
3. Test the pipeline (`nextflow run humann.nf -stub`)