#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -J run_MCScanX
#SBATCH -o run_MCScanX.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# interactive -A uppmax2026-1-8 -t 5:00:00

## runs really fast for five coleopteran species
./MCScanX all_species
