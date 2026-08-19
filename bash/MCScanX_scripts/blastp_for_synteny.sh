#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 5
#SBATCH --mem=30G
#SBATCH -t 18:00:00
#SBATCH -J blastp_for_synteny
#SBATCH -o blastp_for_synteny.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools blast/2.15.0+

## remove the species name prefixes I added for orthofinder
# sed 's/>B_siliquastri__B_siliquastri__B_siliquastri_/>/g' /proj/naiss2023-6-65/Milena/chapter2/protein_data/B_siliquastri.faa > /proj/naiss2023-6-65/Milena/chapter2/protein_data/B_siliquastri_original_header.faa
# sed 's/A_obtectus__A_obtectus__/>/g' /proj/naiss2023-6-65/Milena/chapter2/protein_data/A_obtectus.faa > /proj/naiss2023-6-65/Milena/chapter2/protein_data/A_obtectus_original_header.faa
# sed 's/>C_chinensis__C_chinensis__C_chinensis_/>/g' /proj/naiss2023-6-65/Milena/chapter2/protein_data/C_chinensis.faa > /proj/naiss2023-6-65/Milena/chapter2/protein_data/C_chinensis_original_header.faa
# sed 's/>Cmac_Lome_diverse_/>/g' /proj/naiss2023-6-65/Milena/chapter2/protein_data/C_maculatus_superscaffolded.faa > /proj/naiss2023-6-65/Milena/chapter2/protein_data/C_maculatus_superscaffolded_original_header.faa
# sed 's/>T_castaneum__T_castaneum__T_castaneum_/>/g' /proj/naiss2023-6-65/Milena/chapter2/protein_data/T_castaneum.faa > /proj/naiss2023-6-65/Milena/chapter2/protein_data/T_castaneum_original_header.faa

A_obtectus_proteins=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/A_obtectus.faa
B_siliquastri_proteins=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/B_siliquastri.faa
B_varius_proteins=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/B_varius.faa
C_chinensis_proteins=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/C_chinensis.faa
C_maculatus_proteins=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/C_maculatus.faa
D_carinulata_proteins=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/D_carinulata.faa
D_sublineata_proteins=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/D_sublineata.faa

## --> re-run for new proteinfiles!

## make databases
# for SPECIES1 in $D_sublienata_proteins $D_carinulata_proteins # $A_obtectus_proteins $B_siliquastri_proteins $C_chinensis_proteins $C_maculatus_proteins $T_castaneum_proteins
# do
#     makeblastdb -in $SPECIES1 -dbtype prot
#     echo " ---> done database ${SPECIES1}"
# done
## -->

for SPECIES1 in $B_siliquastri_proteins $B_varius_proteins $C_chinensis_proteins $C_maculatus_proteins $D_carinulata_proteins $D_sublineata_proteins
do  

    SPECIES1_name="${SPECIES1##*/}"
    SPECIES1_name="${SPECIES1_name%.*}"

    for SPECIES2 in $B_siliquastri_proteins $B_varius_proteins $C_chinensis_proteins $C_maculatus_proteins $D_carinulata_proteins $D_sublineata_proteins
    do

        SPECIES2_name="${SPECIES2##*/}"
        SPECIES2_name="${SPECIES2_name%.*}"

        # if [[ "${SPECIES1_name}" == "${SPECIES2_name}" ]]
        # then
        #     continue
        # fi

        OUT_1v2="${SPECIES1_name}_vs_${SPECIES2_name}.blast"
        OUT_2v1="${SPECIES2_name}_vs_${SPECIES1_name}.blast"

        # the documentation says outfmt6 but I think they mean 8
        echo "RUNNING... blastp -query $SPECIES1 -db $SPECIES2 -out $OUT_1v2 -num_threads 5 -num_alignments 5 -evalue 1e-10  -outfmt 6"
        blastp -query $SPECIES1 -db $SPECIES2 -out $OUT_1v2 -num_threads 5 -num_alignments 5 -evalue 1e-10  -outfmt 6
        echo " =========> ${OUT_1v2} done!"

        # # reverse already happens automatically in the nested for loop no need to implement explicitly
        # blastp -query $SPECIES2 -db $SPECIES1 -out $OUT_2v1 -num_threads 5 -num_alignments 5 -evalue 1e-10  -outfmt 6
        # echo " =========> ${OUT_2v1} done!"

    done
done