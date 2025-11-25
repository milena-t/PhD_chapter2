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
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/230/105/GCF_026230105.1_icDioSubl1.1/GCF_026230105.1_icDioSubl1.1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/230/105/GCF_026230105.1_icDioSubl1.1/GCF_026230105.1_icDioSubl1.1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/230/105/GCF_026230105.1_icDioSubl1.1/GCF_026230105.1_icDioSubl1.1_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/230/105/GCF_026230105.1_icDioSubl1.1/GCF_026230105.1_icDioSubl1.1_genomic.fna.gz
gunzip GCF_026230105.1_icDioSubl1.1_protein.faa.gz
gunzip GCF_026230105.1_icDioSubl1.1_genomic.gff.gz
gunzip GCF_026230105.1_icDioSubl1.1_cds_from_genomic.fna.gz
gunzip GCF_026230105.1_icDioSubl1.1_genomic.fna.gz

cd /proj/naiss2023-6-65/Milena/coleoptera_sequences/d_carinulata
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/250/575/GCF_026250575.1_icDioCari1.1/GCF_026250575.1_icDioCari1.1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/250/575/GCF_026250575.1_icDioCari1.1/GCF_026250575.1_icDioCari1.1_genomic.gff.gz 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/250/575/GCF_026250575.1_icDioCari1.1/GCF_026250575.1_icDioCari1.1_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/250/575/GCF_026250575.1_icDioCari1.1/GCF_026250575.1_icDioCari1.1_genomic.fna.gz
gunzip GCF_026250575.1_icDioCari1.1_genomic.gff.gz
gunzip GCF_026250575.1_icDioCari1.1_protein.faa.gz
gunzip GCF_026250575.1_icDioCari1.1_cds_from_genomic.fna.gz
gunzip GCF_026250575.1_icDioCari1.1_genomic.fna.gz