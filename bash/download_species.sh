#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J download_diorhabda
#SBATCH -o download_diorhabda.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

cd /proj/naiss2023-6-65/Milena/coleoptera_sequences/d_sublineata
gunzip GCF_026230105.1_icDioSubl1.1_protein.faa.gz
gunzip GCF_026230105.1_icDioSubl1.1_protein.gff.gz

cd /proj/naiss2023-6-65/Milena/coleoptera_sequences/d_carinulata
gunzip GCF_026250575.1_icDioCari1.1_genomic.gff.gz 
gunzip GCF_026250575.1_icDioCari1.1_protein.faa.gz