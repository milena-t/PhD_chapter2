#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 20
#SBATCH -t 2:00:00
#SBATCH -J orthofinder_chrysomelidae
#SBATCH -o orthofinder_chrysomelidae.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load OrthoFinder/2.5.5-foss-2024a

# the output directory will be created by orthofinder and should NOT already exist!
OUT_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/orthofinder/orthofinder_chrysomelidae

# orthofinder runs with protein sequences here! nucleotide sequences require the -d flag
IN_FASTA_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/orthofinder/proteins_chrysomelidae

# remove it if it does exist, this overwrites all preexisting results
rm -r $OUT_DIR

orthofinder -t 20 -a 20 -f $IN_FASTA_DIR -o $OUT_DIR
