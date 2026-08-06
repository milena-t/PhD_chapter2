#!/bin/sh

ORTHODB_ARTHROPODA=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/annotation_pipeline/annotation_protein_data/OrthoDB_Arthropoda_v11.fa
ASS_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/annotation_pipeline/only_orthodb_annotation

sbatch --job-name="A_obtectus_annotation" --output="A_obtectus_annotation.out" -t 5-00:00:00 \
/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/bash/braker3_singularity_with_RNAseq_in_SNIC_TMP.sh \
A_obtectus ${ASS_DIR}/A_obtectus/assembly_genomic.fna.masked \
$ORTHODB_ARTHROPODA SRR7881623,SRR7881590,SRR7881629,SRR7881599

sbatch --job-name="C_chinensis_annotation" --output="C_chinensis_annotation.out" -t 5-00:00:00 \
/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/bash/braker3_singularity_with_RNAseq_in_SNIC_TMP.sh \
C_chinensis ${ASS_DIR}/C_chinensis/assembly_genomic.fna.masked \
$ORTHODB_ARTHROPODA SRR6167202,SRR6167200,SRR16474572,SRR19159851,SRR19159848,SRR18148876

sbatch --job-name="B_varius_annotation" --output="B_varius_annotation.out" -t 5-00:00:00 \
/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/bash/braker3_singularity_with_RNAseq_in_SNIC_TMP.sh \
B_varius /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/coleoptera_sequences/b_varius/GCA_964204745.1_icBruVari1.hap1.1_genomic.fna \
$ORTHODB_ARTHROPODA ERR14379121