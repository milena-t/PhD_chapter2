#!/bin/bash -l

ASS_DIR=/proj/naiss2023-6-65/Milena/chapter2/assemblies

cd /proj/naiss2023-6-65/Milena/chapter2/pairwise_wga

sbatch --job-name="Diorhabda_WGA" --output="Diorhabda_WGA.out" whole_pipeline.sh "${ASS_DIR}/D_carinulata_assembly.fna" "${ASS_DIR}/D_sublineata_assembly.fna"
sbatch --job-name="Cmac_Aobt_WGA" --output="Cmac_Aobt_WGA.out" whole_pipeline.sh "${ASS_DIR}/C_maculatus_assembly.fna" "${ASS_DIR}/A_obtectus_assembly.fna"
sbatch --job-name="Cmac_Bsil_WGA" --output="Cmac_Bsil_WGA.out" whole_pipeline.sh "${ASS_DIR}/C_maculatus_assembly.fna" "${ASS_DIR}/B_siliquastri_assembly.fna"
sbatch --job-name="Cmac_Dcar_WGA" --output="Cmac_Dcar_WGA.out" whole_pipeline.sh "${ASS_DIR}/C_maculatus_assembly.fna" "${ASS_DIR}/D_carinulata_assembly.fna"
sbatch --job-name="Dcar_Bsil_WGA" --output="Dcar_Bsil_WGA.out" whole_pipeline.sh "${ASS_DIR}/D_carinulata_assembly.fna" "${ASS_DIR}/B_siliquastri_assembly.fna"