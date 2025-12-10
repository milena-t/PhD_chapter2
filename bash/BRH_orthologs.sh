#!/bin/bash -l

BLAST_OUTDIR_UPPMAX=/proj/naiss2023-6-65/Milena/chapter3/all_vs_all_blastp
BLAST_OUTDIR=/Users/miltr339/work/pairwise_blast_chapter_2_3

python3 ../../src/get_blast_BRH.py $BLAST_OUTDIR/A_obtectus_original_header_vs_B_siliquastri_original_header.blast $BLAST_OUTDIR/B_siliquastri_original_header_vs_A_obtectus_original_header.blast