#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 5
#SBATCH --mem=30G
#SBATCH -t 18:00:00
#SBATCH -J blastp_for_synteny
#SBATCH -o blastp_for_synteny.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load BLAST+/2.17.0-gompi-2024a

SPECIES1=$1
SPECIES1_name="${SPECIES1##*/}"
SPECIES1_name="${SPECIES1_name%.*}"

SPECIES2=$2
SPECIES2_name="${SPECIES2##*/}"
SPECIES2_name="${SPECIES2_name%.*}"

OUT_1v2="${SPECIES1_name}_vs_${SPECIES2_name}.blast"

# the documentation says outfmt6 but I think they mean 8
echo "RUNNING... blastp -query $SPECIES1 -db $SPECIES2 -out $OUT_1v2 -num_threads 5 -num_alignments 5 -evalue 1e-10  -outfmt 6"
blastp -query $SPECIES1 -db $SPECIES2 -out $OUT_1v2 -num_threads 5 -num_alignments 5 -evalue 1e-10  -outfmt 6
echo " =========> ${OUT_1v2} done!"
