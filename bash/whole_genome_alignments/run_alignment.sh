#!/bin/bash -l

ASS_DIR=/proj/naiss2023-6-65/Milena/chapter2/assemblies
SCRIPT_DIR=/proj/naiss2023-6-65/Milena/chapter2/PhD_chapter2/bash/whole_genome_alignments

cd /proj/naiss2023-6-65/Milena/chapter2/pairwise_wga

sbatch --job-name="Cmac_Aobt_WGA" --output="Cmac_Aobt_WGA.out" $SCRIPT_DIR/whole_pipeline_rackham.sh "${ASS_DIR}/C_maculatus_superscaffolded_assembly.fna.fna" "${ASS_DIR}/A_obtectus_assembly.fna"
sbatch --job-name="Cmac_Bsil_WGA" --output="Cmac_Bsil_WGA.out" $SCRIPT_DIR/whole_pipeline_rackham.sh "${ASS_DIR}/C_maculatus_superscaffolded_assembly.fna.fna" "${ASS_DIR}/B_siliquastri_assembly.fna"
sbatch --job-name="Diorhabda_WGA" --output="Diorhabda_WGA.out" $SCRIPT_DIR/whole_pipeline_rackham.sh "${ASS_DIR}/D_carinulata_assembly.fna" "${ASS_DIR}/D_sublineata_assembly.fna"
sbatch --job-name="Dcar_Bsil_WGA" --output="Dcar_Bsil_WGA.out" $SCRIPT_DIR/whole_pipeline_rackham.sh "${ASS_DIR}/D_carinulata_assembly.fna" "${ASS_DIR}/B_siliquastri_assembly.fna"
sbatch --job-name="Cmac_Dcar_WGA" --output="Cmac_Dcar_WGA.out" $SCRIPT_DIR/whole_pipeline_rackham.sh "${ASS_DIR}/C_maculatus_superscaffolded_assembly.fna.fna" "${ASS_DIR}/D_carinulata_assembly.fna"