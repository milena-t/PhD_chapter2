#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -n 5
#SBATCH -p core
#SBATCH -t 1:00:00
#SBATCH -J blastp_for_synteny_test_Bsil_Tcas
#SBATCH -o blastp_for_synteny_test_Bsil_Tcas.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools blast/2.15.0+

# blastp -query /proj/naiss2023-6-65/Milena/chapter2/protein_data/B_siliquastri_original_header.faa -db /proj/naiss2023-6-65/Milena/chapter2/protein_data/B_siliquastri_original_header.faa -out /proj/naiss2023-6-65/Milena/chapter2/pairwise_blast_searches_Bsil_Tcas/B_siliquastri_original_header_vs_B_siliquastri_original_header.blast -num_threads 5 -num_alignments 5 -evalue 1e-10  -outfmt 6
# blastp -query /proj/naiss2023-6-65/Milena/chapter2/protein_data/B_siliquastri_original_header.faa -db /proj/naiss2023-6-65/Milena/chapter2/protein_data/T_castaneum_original_header.faa -out /proj/naiss2023-6-65/Milena/chapter2/pairwise_blast_searches_Bsil_Tcas/B_siliquastri_original_header_vs_T_castaneum_original_header.blast -num_threads 5 -num_alignments 5 -evalue 1e-10  -outfmt 6
blastp -query /proj/naiss2023-6-65/Milena/chapter2/protein_data/T_castaneum_original_header.faa -db /proj/naiss2023-6-65/Milena/chapter2/protein_data/B_siliquastri_original_header.faa -out /proj/naiss2023-6-65/Milena/chapter2/pairwise_blast_searches_Bsil_Tcas/T_castaneum_original_header_vs_B_siliquastri_original_header.blast -num_threads 5 -num_alignments 5 -evalue 1e-10  -outfmt 6
# blastp -query /proj/naiss2023-6-65/Milena/chapter2/protein_data/T_castaneum_original_header.faa -db /proj/naiss2023-6-65/Milena/chapter2/protein_data/T_castaneum_original_header.faa -out /proj/naiss2023-6-65/Milena/chapter2/pairwise_blast_searches_Bsil_Tcas/T_castaneum_original_header_vs_T_castaneum_original_header.blast -num_threads 5 -num_alignments 5 -evalue 1e-10  -outfmt 6
