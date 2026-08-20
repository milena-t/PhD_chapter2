#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 5
#SBATCH --mem=30G
#SBATCH -t 18:00:00
#SBATCH -J blastp_for_synteny
#SBATCH -o blastp_for_synteny.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load BLAST+/2.17.0-gompi-2024a

A_obtectus_proteins=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/A_obtectus.faa
B_siliquastri_proteins=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/B_siliquastri.faa
B_varius_proteins=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/B_varius.faa
C_chinensis_proteins=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/C_chinensis.faa
C_maculatus_proteins=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/C_maculatus.faa
D_carinulata_proteins=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/D_carinulata.faa
D_sublineata_proteins=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/D_sublineata.faa

## --> re-run for new proteinfiles!

## make databases
# for SPECIES1 in $A_obtectus_proteins $B_siliquastri_proteins $B_varius_proteins $C_chinensis_proteins $C_maculatus_proteins $D_carinulata_proteins $D_sublineata_proteins
# do
#     makeblastdb -in $SPECIES1 -dbtype prot
#     echo " ---> done database ${SPECIES1}"
# done
## -->

cd /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/blastp_results

for SPECIES1 in $A_obtectus_proteins $B_siliquastri_proteins $B_varius_proteins $C_chinensis_proteins $C_maculatus_proteins $D_carinulata_proteins $D_sublineata_proteins
do  

    for SPECIES2 in $A_obtectus_proteins $B_siliquastri_proteins $B_varius_proteins $C_chinensis_proteins $C_maculatus_proteins $D_carinulata_proteins $D_sublineata_proteins
    do

        SPECIES1_name="${SPECIES1##*/}"
        SPECIES1_name="${SPECIES1_name%.*}"
        SPECIES2_name="${SPECIES2##*/}"
        SPECIES2_name="${SPECIES2_name%.*}"
        OUTFILE="blastp_${SPECIES1_name}_${SPECIES2_name}"
        echo "--------> ${OUTFILE}"
        sbatch -o "${OUTFILE}.out" -J $OUTFILE /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/bash/MCScanX_scripts/blastp_ind_job.sh $SPECIES1 $SPECIES2
        echo ""

    done
done