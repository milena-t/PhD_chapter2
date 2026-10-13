#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH -t 20:00
#SBATCH -J gff_to_fna
#SBATCH -o gff_to_fna.log
#SBATCH --mail-type=ALL


module load gffread/0.12.7-GCCcore-13.3.0 SAMtools/1.22-GCC-13.3.0

ASS_C_chinensis=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/annotation_pipeline/only_orthodb_annotation/C_chinensis/assembly_genomic.fna.masked
ASS_A_obtectus=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/annotation_pipeline/only_orthodb_annotation/A_obtectus/assembly_genomic.fna.masked
ASS_B_varius=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/coleoptera_sequences/b_varius/GCA_964204745.1_icBruVari1.hap1.1_genomic.fna

ANN_C_chinensis=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/braker4_annotation/C_chinensis/results/braker_longest.gff
ANN_A_obtectus=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/braker4_annotation/A_obtectus/results/braker_longest.gff
ANN_B_varius=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/braker4_annotation/B_varius/results/braker_longest.gff

TRA_C_chinensis=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/nucleotide_data/C_chinensis.fna
TRA_A_obtectus=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/nucleotide_data/A_obtectus.fna
TRA_B_varius=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/nucleotide_data/B_varius.fna

samtools faidx $ASS_C_chinensis
samtools faidx $ASS_A_obtectus
samtools faidx $ASS_B_varius


gffread -M -x $TRA_C_chinensis -g $ASS_C_chinensis $ANN_C_chinensis
gffread -M -x $TRA_A_obtectus -g $ASS_A_obtectus $ANN_A_obtectus
gffread -M -x $TRA_B_varius -g $ASS_B_varius $ANN_B_varius