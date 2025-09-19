#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 8:00:00
#SBATCH -J blastp_for_synteny
#SBATCH -o blastp_for_synteny.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools blast/2.15.0+

A_obtectus_proteins=/proj/naiss2023-6-65/Milena/chapter2/protein_data/A_obtectus.faa
B_siliquastri_proteins=/proj/naiss2023-6-65/Milena/chapter2/protein_data/B_siliquastri.faa
C_chinensis_proteins=/proj/naiss2023-6-65/Milena/chapter2/protein_data/C_chinensis.faa
C_maculatus_proteins=/proj/naiss2023-6-65/Milena/chapter2/protein_data/C_maculatus_superscaffolded.faa
T_castaneum_proteins=/proj/naiss2023-6-65/Milena/chapter2/protein_data/T_castaneum.faa

## make databases
for SPECIES1 in $A_obtectus_proteins $B_siliquastri_proteins $C_chinensis_proteins $C_maculatus_proteins $T_castaneum_proteins
do
    makeblastdb -in $SPECIES1 -dbtype prot
    echo " ---> done database ${SPECIES1}"
done

for SPECIES1 in $A_obtectus_proteins $B_siliquastri_proteins $C_chinensis_proteins $C_maculatus_proteins $T_castaneum_proteins
do
    SPECIES1_name="${SPECIES1##*/}"
    SPECIES1_name="${SPECIES1_name%.*}"

    for SPECIES2 in $A_obtectus_proteins $B_siliquastri_proteins $C_chinensis_proteins $C_maculatus_proteins $T_castaneum_proteins
    do

        SPECIES2_name="${SPECIES2##*/}"
        SPECIES2_name="${SPECIES2_name%.*}"

        OUT_1v2="${SPECIES1_name}_vs_${SPECIES2_name}.blast"
        OUT_2v1="${SPECIES2_name}_vs_${SPECIES1_name}.blast"

        blastp -query $SPECIES1 -db $SPECIES2 -out $OUT_1v2 -num_threads 5 -max_target_seqs 5 -evalue 1e-10  -outfmt 8
        echo " =========> ${OUT_1v2} done!"

        blastp -query $SPECIES2 -db $SPECIES1 -out $OUT_2v1 -num_threads 5 -max_target_seqs 5 -evalue 1e-10  -outfmt 8
        echo " =========> ${OUT_2v1} done!"

    done
done