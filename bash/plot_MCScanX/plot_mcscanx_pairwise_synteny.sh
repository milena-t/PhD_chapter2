#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -n 1
#SBATCH -t 10:00
#SBATCH -J plot_MCScanX
#SBATCH -o plot_MCScanX.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# Bioinfo tools not needed on pelle

cd /proj/naiss2023-6-65/Milena/chapter2/MCScanX-1_0_0/downstream_analyses

java dual_synteny_plotter \
    -g /proj/naiss2023-6-65/Milena/chapter2/MCScanX_results/all_species.gff \
    -s /proj/naiss2023-6-65/Milena/chapter2/MCScanX_results/all_species.collinearity \
    -c plot_ctl_Tcas_Bsil.ctl \
    -o /proj/naiss2023-6-65/Milena/chapter2/MCScanX_results/plots/plot_ctl_Tcas_Bsil.png

