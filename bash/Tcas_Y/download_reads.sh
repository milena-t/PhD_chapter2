#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -n 1
#SBATCH -p core
#SBATCH -t 1:00:00
#SBATCH -J SRR23732241_download_sequences
#SBATCH -o SRR23732241_download_sequences.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# Bioinfo tools not needed on pelle
# module load bioinfo-tools EDirect



## use with SRR IDS from PhD_chapter2/bash/Tcas_Y/SRR_list.sh
curl -O https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=$1