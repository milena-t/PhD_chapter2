#!/bin/bash -l

ASS_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/assemblies
SCRIPT_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/bash/whole_genome_alignments

cd /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/pairwise_wga/systematic_pairs

A_obtectus=A_obtectus_assembly.fna
B_siliquastri=B_siliquastri.masked.fna
B_varius=B_varius_assembly.fna
C_maculatus=C_maculatus.masked.fna
C_chinensis=C_chinensis.masked.fna
D_carinulata=D_carinulata_assembly.fna
D_sublineata=D_sublineata_assembly.fna

sbatch --job-name="A_obtectus_B_siliquastri_WGA" --output="A_obtectus_B_siliquastri_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${A_obtectus}" "${ASS_DIR}/${B_siliquastri}"
sbatch --job-name="A_obtectus_B_varius_WGA" --output="A_obtectus_B_varius_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${A_obtectus}" "${ASS_DIR}/${B_varius}"
sbatch --job-name="A_obtectus_C_maculatus_WGA" --output="A_obtectus_C_maculatus_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${A_obtectus}" "${ASS_DIR}/${C_maculatus}"
sbatch --job-name="A_obtectus_C_chinensis_WGA" --output="A_obtectus_C_chinensis_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${A_obtectus}" "${ASS_DIR}/${C_chinensis}"
sbatch --job-name="A_obtectus_D_carinulata_WGA" --output="A_obtectus_D_carinulata_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${A_obtectus}" "${ASS_DIR}/${D_carinulata}"
sbatch --job-name="A_obtectus_D_sublineata_WGA" --output="A_obtectus_D_sublineata_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${A_obtectus}" "${ASS_DIR}/${D_sublineata}"

sbatch --job-name="B_siliquastri_B_varius_WGA" --output="B_siliquastri_B_varius_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${B_siliquastri}" "${ASS_DIR}/${B_varius}"
sbatch --job-name="B_siliquastri_C_maculatus_WGA" --output="B_siliquastri_C_maculatus_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${B_siliquastri}" "${ASS_DIR}/${C_maculatus}"
sbatch --job-name="B_siliquastri_C_chinensis_WGA" --output="B_siliquastri_C_chinensis_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${B_siliquastri}" "${ASS_DIR}/${C_chinensis}"
sbatch --job-name="B_siliquastri_D_carinulata_WGA" --output="B_siliquastri_D_carinulata_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${B_siliquastri}" "${ASS_DIR}/${D_carinulata}"
sbatch --job-name="B_siliquastri_D_sublineata_WGA" --output="B_siliquastri_D_sublineata_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${B_siliquastri}" "${ASS_DIR}/${D_sublineata}"

sbatch --job-name="B_varius_C_maculatus_WGA" --output="B_varius_C_maculatus_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${B_varius}" "${ASS_DIR}/${C_maculatus}"
sbatch --job-name="B_varius_C_chinensis_WGA" --output="B_varius_C_chinensis_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${B_varius}" "${ASS_DIR}/${C_chinensis}"
sbatch --job-name="B_varius_D_carinulata_WGA" --output="B_varius_D_carinulata_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${B_varius}" "${ASS_DIR}/${D_carinulata}"
sbatch --job-name="B_varius_D_sublineata_WGA" --output="B_varius_D_sublineata_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${B_varius}" "${ASS_DIR}/${D_sublineata}"

sbatch --job-name="C_maculatus_C_chinensis_WGA" --output="C_maculatus_C_chinensis_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${C_maculatus}" "${ASS_DIR}/${C_chinensis}"
sbatch --job-name="C_maculatus_D_carinulata_WGA" --output="C_maculatus_D_carinulata_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${C_maculatus}" "${ASS_DIR}/${D_carinulata}"
sbatch --job-name="C_maculatus_D_sublineata_WGA" --output="C_maculatus_D_sublineata_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${C_maculatus}" "${ASS_DIR}/${D_sublineata}"

sbatch --job-name="C_chinensis_D_carinulata_WGA" --output="C_chinensis_D_carinulata_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${C_chinensis}" "${ASS_DIR}/${D_carinulata}"
sbatch --job-name="C_chinensis_D_sublineata_WGA" --output="C_chinensis_D_sublineata_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${C_chinensis}" "${ASS_DIR}/${D_sublineata}"

sbatch --job-name="D_carinulata_D_sublineata_WGA" --output="D_carinulata_D_sublineata_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${D_carinulata}" "${ASS_DIR}/${D_sublineata}"

# sbatch --job-name="Cmac_Aobt_WGA" --output="Cmac_Aobt_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${C_maculatus}" "${ASS_DIR}/${A_obtectus}"
# sbatch --job-name="Cmac_Bsil_WGA" --output="Cmac_Bsil_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${C_maculatus}" "${ASS_DIR}/${B_siliquastri}"
# sbatch --job-name="Cmac_Dcar_WGA" --output="Cmac_Dcar_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${C_maculatus}" "${ASS_DIR}/${D_carinulata}"
# sbatch --job-name="Diorhabda_WGA" --output="Diorhabda_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${D_carinulata}" "${ASS_DIR}/${D_sublineata}"
# sbatch --job-name="Dcar_Bsil_WGA" --output="Dcar_Bsil_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${D_carinulata}" "${ASS_DIR}/${B_siliquastri}"
# sbatch --job-name="Bruchidius_WGA" --output="Bruchidius_WGA.out" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${B_varius}" "${ASS_DIR}/${B_siliquastri}"