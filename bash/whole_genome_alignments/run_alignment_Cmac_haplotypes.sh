#!/bin/bash -l


CMAC_Y_s=/proj/coleoptera-genomics-2025/snic2021-6-30/Philipp/HiFi/Whole_genome_alinment/Alinment_pt_036/pt_036_001_hifiasm_primary.fasta
CMAC_Y_l=/proj/coleoptera-genomics-2025/snic2021-6-30/Philipp/HiFi/Whole_genome_alinment/Alinment_pt_036/pt_036_002_hifiasm_primary.fasta


for SPECIES1 in $CMAC_Y_s $CMAC_Y_l
do  

    SPECIES1_name="${SPECIES1##*/}"
    SPECIES1_name="${SPECIES1_name%.masked.fna*}"

    for SPECIES2 in $CMAC_Y_s $CMAC_Y_l
    do

        SPECIES2_name="${SPECIES2##*/}"
        SPECIES2_name="${SPECIES2_name%.masked.fna*}"

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