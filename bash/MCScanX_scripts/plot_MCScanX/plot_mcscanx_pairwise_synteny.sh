#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH -t 10:00
#SBATCH -J plot_MCScanX
#SBATCH -o plot_MCScanX.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# Bioinfo tools not needed on pelle

FILES_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/

cd /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/MCScanX-1_0_0/downstream_analyses


java dual_synteny_plotter \
    -g ${FILES_DIR}MCScanX_results/all_species.gff \
    -s ${FILES_DIR}MCScanX_results/all_species.collinearity \
    -c ${FILES_DIR}PhD_chapter2/bash/MCScanX_scripts/plot_MCScanX/plot_ctl_Cmac_Bsil.ctl \
    -o ${FILES_DIR}MCScanX_results/plots/plot_ctl_Cmac_Bsil.png

echo "--> done Cmac Bsil"

java dual_synteny_plotter \
    -g ${FILES_DIR}MCScanX_results/all_species.gff \
    -s ${FILES_DIR}MCScanX_results/all_species.collinearity \
    -c ${FILES_DIR}PhD_chapter2/bash/MCScanX_scripts/plot_MCScanX/plot_ctl_Aobt_Bsil.ctl \
    -o ${FILES_DIR}MCScanX_results/plots/plot_ctl_Aobt_Bsil.png

echo "--> done Aobt Bsil"

java dual_synteny_plotter \
    -g ${FILES_DIR}MCScanX_results/all_species.gff \
    -s ${FILES_DIR}MCScanX_results/all_species.collinearity \
    -c ${FILES_DIR}PhD_chapter2/bash/MCScanX_scripts/plot_MCScanX/plot_ctl_Aobt_Bvar.ctl \
    -o ${FILES_DIR}MCScanX_results/plots/plot_ctl_Aobt_Bvar.png

echo "--> done Aobt Bvar"

java dual_synteny_plotter \
    -g ${FILES_DIR}MCScanX_results/all_species.gff \
    -s ${FILES_DIR}MCScanX_results/all_species.collinearity \
    -c ${FILES_DIR}PhD_chapter2/bash/MCScanX_scripts/plot_MCScanX/plot_ctl_Dsub_Dcar.ctl \
    -o ${FILES_DIR}MCScanX_results/plots/plot_ctl_Dsub_Dcar.png

echo "--> done Dsub Dcar"

java dual_synteny_plotter \
    -g ${FILES_DIR}MCScanX_results/all_species.gff \
    -s ${FILES_DIR}MCScanX_results/all_species.collinearity \
    -c ${FILES_DIR}PhD_chapter2/bash/MCScanX_scripts/plot_MCScanX/plot_ctl_Bvar_Bsil.ctl \
    -o ${FILES_DIR}MCScanX_results/plots/plot_ctl_Bvar_Bsil.png

echo "--> done Bvar Bsil"

java dual_synteny_plotter \
    -g ${FILES_DIR}MCScanX_results/all_species.gff \
    -s ${FILES_DIR}MCScanX_results/all_species.collinearity \
    -c ${FILES_DIR}PhD_chapter2/bash/MCScanX_scripts/plot_MCScanX/plot_ctl_Cmac_Bvar.ctl \
    -o ${FILES_DIR}MCScanX_results/plots/plot_ctl_Cmac_Bvar.png

echo "--> done Cmac Bvar"