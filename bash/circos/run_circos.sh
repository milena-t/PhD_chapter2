#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -J run_circos
#SBATCH -o run_circos.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# interactive -A uppmax2026-1-8 -t 5:00:00
module load Circos/0.69-10-GCCcore-13.3.0

python3 /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/src/circos/loop_circos_species.py


