#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 2:00:00
#SBATCH -J orthofinder_TE_orthoDB
#SBATCH -o orthofinder_TE_filtered_orthoDB.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools OrthoFinder/2.5.2

# the output directory will be created by orthofinder and should NOT already exist!
OUT_DIR=/proj/naiss2023-6-65/Milena/chapter2/orthofinder/orthofinder_Cmac_Kaufmann2023

# orthofinder runs with protein sequences here! nucleotide sequences require the -d flag
IN_FASTA_DIR=/proj/naiss2023-6-65/Milena/chapter2/protein_data/orthofinder_dir

# remove it if it does exist, this overwrites all preexisting results
rm -r $OUT_DIR

orthofinder.py -t 20 -a 20 -f $IN_FASTA_DIR -o $OUT_DIR
