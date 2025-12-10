#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 1:00:00
#SBATCH -J blastp_for_synteny
#SBATCH -o blastp_for_synteny.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools blast/2.15.0+

SPECIES1=/proj/naiss2023-6-65/Milena/chapter2/protein_data/T_castaneum_original_header.faa
SPECIES2=/proj/naiss2023-6-65/Milena/chapter2/protein_data/C_maculatus_superscaffolded_original_header.faa

SPECIES1_name="${SPECIES1##*/}"
SPECIES1_name="${SPECIES1_name%.*}"
SPECIES2_name="${SPECIES2##*/}"
SPECIES2_name="${SPECIES2_name%.*}"

OUT_1v2="${SPECIES1_name}_vs_${SPECIES2_name}.blast"

blastp -query $SPECIES1 -db $SPECIES2 -out $OUT_1v2 -num_threads 5 -max_target_seqs 5 -evalue 1e-10  -outfmt 6
echo " =========> ${OUT_1v2} done!"