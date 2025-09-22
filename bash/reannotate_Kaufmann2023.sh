#!/bin/sh

ORTHODB_ARTHROPODA=/proj/naiss2023-6-65/Milena/annotation_pipeline/annotation_protein_data/OrthoDB_Arthropoda_v11.fa

LOME_DIR=/proj/naiss2023-6-65/Milena/annotation_pipeline/Kaufmann2023_updated_RNAseq_annotation/Cmac_not_superscaffolded
cd $LOME_DIR
sbatch --job-name="Kaufmann2023_reannotate" --output="Kaufmann2023_reannotate.out" -t 5-00:00:00 \
run_braker.sh Cmac_Kaufmann2023 ${LOME_DIR}/assembly_genomic.fna.masked \
$ORTHODB_ARTHROPODA ERR12383274,ERR12383248,ERR12383255,ERR12383289,ERR12383294,ERR12383314,ERR12383271,ERR12383290,ERR12383279,ERR12383266