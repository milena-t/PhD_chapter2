#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -n 5
#SBATCH -t 1:00:00
#SBATCH -J read_transcripts_from_gff
#SBATCH -o read_transcripts_from_gff.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# use gffread to extract the protein coding sequences
# -M :  cluster the input transcripts into loci, discarding "duplicated" transcripts (those with the same exact introns and fully contained or equal boundaries)
# -x :  write a FASTA file with spliced CDS for each GFF transcript

module load gffread/0.12.7-GCCcore-13.3.0 SAMtools/1.22-GCC-13.3.0 AGAT/1.6.1-GCCcore-13.3.0

ASS_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/assemblies
MINI_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/paralogs_tblastn/miniprot

for SPECIES in A_obtectus B_siliquastri B_varius C_chinensis C_maculatus D_carinulata D_sublineata 

do 
    echo " *** -${SPECIES}- *** "

    ANNOT_GFF_RAW="${MINI_DIR}/${SPECIES}.masked_align.gff"
    ANNOT_GFF="${MINI_DIR}/${SPECIES}.masked_align_noPAF.gff"
    ASSEMBLY="${ASS_DIR}/${SPECIES}.masked.fna" 

    echo "---------------------------------------------"
    echo " * ${ANNOT_GFF_RAW}"
    echo " * ${ASSEMBLY}"
    echo "---------------------------------------------"

    ANNOT_TRANSCRIPTS=${ANNOT_GFF_RAW}_transcripts.fna
    ANNOT_PROTEINS=${ANNOT_GFF_RAW}_proteins.faa
    TRANSEQ_PATH=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/software_install/emboss/EMBOSS-6.6.0/emboss/transeq

    echo $(pwd)
    echo $(ls -lh $ASSEMBLY)

    ## fix IDs and cds reading frame with AGAT
    # ANNOT_GFF_ID=${ANNOT_GFF_RAW}_AGAT_ID.gff
    # agat_sp_manage_IDs.pl --gff $ANNOT_GFF_RAW -o $ANNOT_GFF_ID
    # ANNOT_GFF=${ANNOT_GFF_ID}_CDS.gff
    # agat_sp_fix_cds_phases.pl --gff $ANNOT_GFF_ID --fasta $ASSEMBLY -o $ANNOT_GFF
    # rm $ANNOT_GFF_ID

    # index assemblies (greatly decreases computing time, and won't work for the more fragmented callosobruchus assemblies otherwise)
    # samtools faidx $ASSEMBLY

    ## remove PAF lines
    grep -v '^##PAF' $ANNOT_GFF_RAW > $ANNOT_GFF

    # extract transcript sequences
    echo "gffread $ANNOT_GFF -M -x $ANNOT_TRANSCRIPTS -g $ASSEMBLY"
    gffread $ANNOT_GFF -M -x $ANNOT_TRANSCRIPTS -g $ASSEMBLY

    # gffread $ANNOT_GFF -g $ASSEMBLY -y $ANNOT_PROTEINS

    # change fasta headers to include species names
    # sed -i "s/>/>${SPECIES_NAME}_/g" $ANNOT_TRANSCRIPTS
    # translate transcript sequences
    $TRANSEQ_PATH -sequence $ANNOT_TRANSCRIPTS -outseq $ANNOT_PROTEINS
    ls -lh $ANNOT_TRANSCRIPTS
    echo "###########################################"

done