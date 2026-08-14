#!/bin/bash -l

#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -t 2-00:00:00
#SBATCH -J braker4_snakemake_Bvar_%j
#SBATCH -o braker4_snakemake_Bvar_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se


module load snakemake/8.27.0-foss-2024a
cd /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/braker4_annotation

# echo " * Augustus config does not exist, create it and change write permissions"
# module load AUGUSTUS/3.5.0-gfbf-2024a # so that the source command works
# cp -dR --preserve=mode,timestamps --no-preserve=ownership $AUGUSTUS_CONFIG_PATH AUGUSTUS_config
# chmod -R +w AUGUSTUS_config
# module unload AUGUSTUS/3.5.0-gfbf-2024a # same as above, some weird shit with conflicting perl versions
# export AUGUSTUS_CONFIG_PATH=$PWD/AUGUSTUS_config

snakemake \
    --keep-going \
    --snakefile /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/software_install/BRAKER4/Snakefile \
    --executor slurm \
    --default-resources slurm_partition=pelle mem_mb=120000 slurm_account=uppmax2026-1-8\
    --cores 20 \
    --jobs 50 \
    --use-singularity \
    --singularity-prefix .singularity_cache \
    --singularity-args "-B /home -B /scratch -B /gorilla -B /proj" \
    --latency-wait 120 \
    --restart-times 3