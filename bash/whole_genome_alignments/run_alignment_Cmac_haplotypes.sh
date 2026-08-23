#!/bin/bash -l


ASS_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/pairwise_wga/Cmac_populations/assemblies
SCRIPT_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/bash/whole_genome_alignments

CMAC_Y_s=Cmac_china.fasta
CMAC_Y_l=Lome_Yl.fasta
CMAC_CHINA=Lome_Ys.fasta

cd /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/pairwise_wga/Cmac_populations

for SPECIES1 in $CMAC_Y_s $CMAC_Y_l $CMAC_CHINA
do  

    SPECIES1_name="${SPECIES1##*/}"
    SPECIES1_name="${SPECIES1_name%.fasta*}"

    for SPECIES2 in $CMAC_Y_s $CMAC_Y_l $CMAC_CHINA
    do

        SPECIES2_name="${SPECIES2##*/}"
        SPECIES2_name="${SPECIES2_name%.fasta*}"

        if [[ "${SPECIES1_name}" == "${SPECIES2_name}" ]]
        then
            continue
        fi

        logfile="${SPECIES1_name}_${SPECIES2_name}_WGA.out"
        jobname="${SPECIES1_name}_${SPECIES2_name}_WGA"
        echo "species1 : ${SPECIES1_name}"
        echo "species2 : ${SPECIES2_name}"
        echo "sbatch --job-name=${jobname} --output=${logfile} $SCRIPT_DIR/whole_pipeline.sh ${ASS_DIR}/${SPECIES1} ${ASS_DIR}/${SPECIES2}"

        sbatch --job-name="${jobname}" --output="${logfile}" $SCRIPT_DIR/whole_pipeline.sh "${ASS_DIR}/${SPECIES1}" "${ASS_DIR}/${SPECIES2}"
        echo ""

    done
done