#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 16
#SBATCH --mem=30G
#SBATCH -t 18:00:00
#SBATCH -J miniprot_for_paralogs
#SBATCH -o miniprot_for_paralogs.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load miniprot/0.18-GCCcore-13.3.0

PROTEINS=$1
ASSEMBLY=$2

SPECIES="${ASSEMBLY##*/}"
SPECIES="${SPECIES%.*}"

echo "...running... minisplice predict -t16 -c vi2-7k.kan.cali vi2-7k.kan ${ASSEMBLY} > ${SPECIES}.tsv"
minisplice predict -t16 -c vi2-7k.kan.cali vi2-7k.kan ${ASSEMBLY} > ${SPECIES}.tsv
echo " =========> minisplice ${OUTFILE} done!"

echo "...running... miniprot -I -G 500000 -u -t16 --gff -j2 --spsc=${SPECIES}.tsv ${ASSEMBLY} ${PROTEINS} > ${SPECIES}_align.gff"
miniprot -I -G 500000 -u -t16 --gff -j2 --spsc=${SPECIES}.tsv ${ASSEMBLY} ${PROTEINS} > "${SPECIES}_align.gff"
echo " =========> miniprot ${OUTFILE} done!"
