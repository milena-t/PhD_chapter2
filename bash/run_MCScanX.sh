#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J run_MCScanX
#SBATCH -o run_MCScanX.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# Bioinfo tools not needed on pelle
# module load bioinfo-tools 

## runs really fast for five coleopteran species
/proj/naiss2023-6-65/Milena/chapter2/MCScanX-1.0.0/MCScanX /proj/naiss2023-6-65/Milena/chapter2/MCScanX_results