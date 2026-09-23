#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 5
#SBATCH --mem=30G
#SBATCH -t 18:00:00
#SBATCH -J tblastn_for_paralogs
#SBATCH -o tblastn_for_paralogs.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load BLAST+/2.17.0-gompi-2024a

PROTEINS=$1
ASSEMBLY=$2
OUTFILE=$3

# the documentation says outfmt6 but I think they mean 8
echo "RUNNING... tblastn -query $PROTEINS -db $ASSEMBLU -out $OUTFILE -num_threads 5 -max_intron_length 3000 -evalue 1e-3 -outfmt 6"
tblastn -query $PROTEINS -db $ASSEMBLU -out $OUTFILE -num_threads 5 -max_intron_length 3000 -evalue 1e-3 -outfmt 6
echo " =========> ${OUTFILE} done!"
