#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J blastp_for_synteny
#SBATCH -o blastp_for_synteny.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools 

./proj/naiss2023-6-65/Milena/chapter2/MCScanX-1.0.0/MCScanX /proj/naiss2023-6-65/Milena/chapter2/MCScanX_results