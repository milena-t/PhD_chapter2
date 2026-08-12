#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 8
#SBATCH -t 1:00:00
#SBATCH -J download_sequences_Cchi
#SBATCH -o download_sequences_Cchi.log
#SBATCH --mail-type=ALL

module load pigz/2.8-GCCcore-13.3.0

# Define directories
BASE=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/braker4_annotation/rnaseq_downloads
SRA=$BASE/B_varius/sra
FASTQ=$BASE/B_varius/fastq
SRAPATH=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/software_install/sra_tools/sratoolkit.3.3.0-ubuntu64/bin/

mkdir -p $SRA $FASTQ


SRR_IDS=(
    SRR6167202
    SRR6167200
    SRR16474572
    SRR19159851
    SRR19159848
    SRR18148876
)

for acc in "${SRR_IDS[@]}"; do
    echo "================== ${acc} =================="
    ${SRAPATH}prefetch $acc --output-directory $SRA # Download .sra file
    ${SRAPATH}fasterq-dump --split-files --threads 4 "$acc"
    gzip "${acc}"*.fastq
    # rm -rf "$acc"
done

