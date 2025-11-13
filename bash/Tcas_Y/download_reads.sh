#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -n 1
#SBATCH -p core
#SBATCH -t 1:00:00
#SBATCH -J SRR23732240_download_sequences
#SBATCH -o SRR23732240_download_sequences.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# Bioinfo tools not needed on pelle
module load bioinfo-tools EDirect

## SRR IDS
curl -O https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR23732240
# curl -O https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR23732241
# curl -O https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR23732242
# curl -O https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR23732243
# curl -O https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR23732244
# curl -O https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR23732245