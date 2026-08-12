#!/bin/bash -l

module load snakemake/8.27.0-foss-2024a
cd /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/assemblies

snakemake \
    --snakefile /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/software_install/BRAKER4/Snakefile \
    --executor slurm \
    --default-resources slurm_partition=pelle mem_mb=120000 \
    --cores 20 \
    --jobs 50 \
    --use-singularity \
    --singularity-prefix .singularity_cache \
    --singularity-args "-B /home -B /scratch" \
    --latency-wait 120 \
    --restart-times 3