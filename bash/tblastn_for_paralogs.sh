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

ASS_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/assemblies
A_obtectus=${ASS_DIR}/A_obtectus.masked.fna
B_siliquastri=${ASS_DIR}/B_siliquastri.masked.fna
B_varius=${ASS_DIR}/B_varius.masked.fna
C_chinensis=${ASS_DIR}/C_chinensis.masked.fna
C_maculatus=${ASS_DIR}/C_maculatus.masked.fna
D_carinulata=${ASS_DIR}/D_carinulata.masked.fna
D_sublineata=${ASS_DIR}/D_sublineata.masked.fna

PROT_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data
A_obtectus_proteins=${PROT_DIR}/A_obtectus.faa
B_siliquastri_proteins=${PROT_DIR}/B_siliquastri.faa
B_varius_proteins=${PROT_DIR}/B_varius.faa
C_chinensis_proteins=${PROT_DIR}/C_chinensis.faa
C_maculatus_proteins=${PROT_DIR}/C_maculatus.faa
D_carinulata_proteins=${PROT_DIR}/D_carinulata.faa
D_sublineata_proteins=${PROT_DIR}/D_sublineata.faa

OUTDIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/paralogs_tblastn

## make databases
cd $ASS_DIR

for SPECIES in "A_obtectus" "B_siliquastri" "B_varius" "C_chinensis" "C_maculatus" "D_carinulata" "D_sublineata"
do
    ASSEMBLY="${ASS_DIR}/${SPECIES}.masked.fna"
    makeblastdb -in $SPECIES1 -dbtype nucl
    echo " ---> done database ${ASSEMBLY}"
done
## -->

cd $OUTDIR

for SPECIES in "A_obtectus" "B_siliquastri" "B_varius" "C_chinensis" "C_maculatus" "D_carinulata" "D_sublineata"
do  
    PROTEINS="${PROT_DIR}/${SPECIES}.faa"
    ASSEMBLY="${ASS_DIR}/${SPECIES}.masked.fna"
    OUTFILE="tblastn_${SPECIES}.out"
    
    echo "--------> ${OUTFILE}"
    sbatch -o "${OUTFILE}" -J $OUTFILE /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/bash/tblastn_ind_job.sh $PROTEINS $ASSEMBLY $OUTFILE
    echo ""

done