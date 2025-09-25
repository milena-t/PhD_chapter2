#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -n 1
#SBATCH -p core
#SBATCH -t 1:00:00
#SBATCH -J run_GENESPACE
#SBATCH -o run_GENESPACE.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# Bioinfo tools not needed on pelle
module load bioinfo-tools OrthoFinder/2.5.2 R/4.2.1

Rscript ../src/genespace.R