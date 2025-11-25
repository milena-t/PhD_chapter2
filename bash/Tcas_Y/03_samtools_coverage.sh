#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J samtools_coverages
#SBATCH -o samtools_coverages.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools samtools/1.20

MAP_DIR=/proj/naiss2023-6-65/Milena/chapter2/Tribolium_poolseq/mapped_data

samtools coverage -b /proj/naiss2023-6-65/Milena/chapter2/PhD_chapter2/bash/Tcas_Y/samtools_list.sh -D -o /proj/naiss2023-6-65/Milena/chapter2/Tribolium_poolseq/contig_coverage.txt