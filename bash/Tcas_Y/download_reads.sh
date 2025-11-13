#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -n 8
#SBATCH -p core
#SBATCH -t 10:00:00
#SBATCH -J SRR23732241_download_sequences
#SBATCH -o SRR23732241_download_sequences.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools
module load sratools/3.0.7
module load pigz

# Define directories
BASE=/proj/naiss2023-6-65/Milena/chapter2/Tribolium_poolseq
SRA=$BASE/raw_data/sra
FASTQ=$BASE/raw_data/fastq

mkdir -p $SRA $FASTQ

prefetch $1 --output-directory $SRA # Download .sra file
fasterq-dump $SRA/$run/$run.sra -O $FASTQ --split-files --threads 8 # Convert to FASTQ (paired-end, gzipped)

pigz -p 8 $FASTQ/*.fastq # Compress