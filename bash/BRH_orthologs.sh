#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -J list_sex_chr_BRHs
#SBATCH -o list_sex_chr_BRHs.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

BLAST_OUTDIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/blastp_results
ANN_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/native_annotations
SCRIPT_PATH=/Users/miltr339/work/PhD_code/PhD_chapter2/src


for SPECIES1 in A_obtectus B_siliquastri B_varius C_chinensis C_maculatus D_carinulata D_sublineata 
do 
    for SPECIES2 in A_obtectus B_siliquastri B_varius C_chinensis C_maculatus D_carinulata D_sublineata 
    do 

        [[ "$SPECIES1" == "$SPECIES2" ]] && continue # skip self-blast

        echo "------------------ ${SPECIES1} vs. ${SPECIES2} ------------------"
        
        python3 $SCRIPT_PATH/get_blast_BRH.py \
        --blast1 ${BLAST_OUTDIR}/${SPECIES1}_vs_${SPECIES2}.blast \
        --blast2 ${BLAST_OUTDIR}/${SPECIES2}_vs_${SPECIES1}.blast \
        --annotation1 $ANN_DIR/${SPECIES1}.gff \
        --annotation2 $ANN_DIR/${SPECIES2}.gff 

    done
done
