#!/bin/bash -l

ASS_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/assemblies
SCRIPT_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/bash/whole_genome_alignments

cd /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/pairwise_wga

sbatch --job-name="Cmac_Aobt_WGA" --output="Cmac_Aobt_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/C_maculatus_superscaffolded_assembly.fna" "${ASS_DIR}/A_obtectus_assembly.fna"
sbatch --job-name="Cmac_Bsil_WGA" --output="Cmac_Bsil_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/C_maculatus_superscaffolded_assembly.fna" "${ASS_DIR}/B_siliquastri_assembly.fna"
sbatch --job-name="Cmac_Dcar_WGA" --output="Cmac_Dcar_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/C_maculatus_superscaffolded_assembly.fna" "${ASS_DIR}/D_carinulata_assembly.fna"
sbatch --job-name="Diorhabda_WGA" --output="Diorhabda_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/d_carinulata_assembly.fna" "${ASS_DIR}/d_sublineata_assembly.fna"
sbatch --job-name="Dcar_Bsil_WGA" --output="Dcar_Bsil_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/d_carinulata_assembly.fna" "${ASS_DIR}/B_siliquastri_assembly.fna"
sbatch --job-name="Bruchidius_WGA" --output="Bruchidius_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/B_varius_assembly.fna" "${ASS_DIR}/B_siliquastri_assembly.fna"