#!/bin/bash -l

BLAST_OUTDIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/blastp_results
ANN_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/native_annotations
SCRIPT_PATH=/Users/miltr339/work/PhD_code/PhD_chapter2/src


# python3 $SCRIPT_PATH/get_blast_BRH.py \
# --blast1 $BLAST_OUTDIR/T_castaneum_original_header_vs_T_freemani_original_header.blast \
# --blast2 $BLAST_OUTDIR/T_freemani_original_header_vs_T_castaneum_original_header.blast \
# --annotation1 $ANN_DIR/T_castaneum.gff \
# --annotation2 $ANN_DIR/T_freemani.gff \
# --X_contigs1 NC_087403.1 \
# --X_contigs2 ENA_CM039461_CM039461.1_Tribolium_freemani_isolate_YK \
