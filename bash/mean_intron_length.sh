#!/bin/bash -l

ASS_DIR=/Users/miltr339/work/chapter2/native_annotations

## filter for only introns with grep and then calculate length by abs($4-$5) with awk


for SPECIES in  "A_obtectus" "B_siliquastri" "B_varius" "C_chinensis" "C_maculatus" "D_carinulata" "D_sublineata" 
do
    ASSEMBLY="${ASS_DIR}/${SPECIES}.gff"
    echo "-------------- ${SPECIES} --------------"
    grep "intron" $ASSEMBLY | awk -F'\t' '{ diff = $4 - $5; sum += (diff < 0 ? -diff : diff) } END { print (NR > 0 ? sum / NR : 0) }'
done 