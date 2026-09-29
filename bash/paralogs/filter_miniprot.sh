#!/bin/bash -l

MINI_DIR=/Users/miltr339/work/chapter2/miniprot_annot
ANN_DIR=/Users/miltr339/work/chapter2/native_annotations

#### FILTER MINIPROT OUTPUT TO GET TRUE AMPLICON GROUPS
# 1. filter for only mapped regions that aren't overlapping with annotated genes, true non-annotated paralogs/pseudogenes
# 2. filter regions that have multiple miniprot aln hits to retain only the best one 
#    (only one query can align to any place in the genome, so each region can only be an amplicon to one source)

cd $MINI_DIR

for SPECIES in  "A_obtectus" "B_siliquastri" "B_varius" "C_chinensis" "C_maculatus" "D_carinulata" "D_sublineata" 
do
    MINIPROT="${MINI_DIR}/${SPECIES}.masked_align.gff"
    ANNOTATION="${ANN_DIR}/${SPECIES}.gff"
    echo "-------------- ${SPECIES} --------------"
    ## get geneIDs that should be removed

    # only mRNA features from miniprot
    grep -v '^#' "${MINIPROT}" | awk '$3=="mRNA"' > mp_mRNA.gff
    # only transcript features from annotation 
    if [ "${SPECIES}" = "D_carinulata" ] || [ "${SPECIES}" = "D_sublineata" ] ;then
        awk '$3=="mRNA"' "${ANNOTATION}" > genes.gff
    else
        awk '$3=="transcript"' "${ANNOTATION}" > genes.gff
    fi
    
    echo "annotated genes:"
    wc -l genes.gff
    echo "miniprot hits before filtering:"
    wc -l mp_mRNA.gff

    ## intersect by position, keep IDs where a mapped miniprot region intersects with annotated genes (this also excludes all self-hits)
    ANNOT_INT_LIST="${SPECIES}_annotation_self_intersecting.txt"
    # only self-strand overlaps with -s
    bedtools intersect -wa -wb -s -a mp_mRNA.gff -b genes.gff | awk -F'\t' '{
        match($9,/(^|;)ID=[^;]+/); id=substr($9,RSTART,RLENGTH); sub(/^;?ID=/,"",id);
        print id
    }' | sort -u > "${ANNOT_INT_LIST}"

    ## remove all IDs in the list generated above
    # this does not include the #PAF lines because gffread does not like them
    NOCROSS=${SPECIES}_miniprot_no_cross_hits.gff
    grep -v '^##PAF' "${MINIPROT}" | gffread - -F --nids "${ANNOT_INT_LIST}" -o "${NOCROSS}"

    echo "   alignent IDs that intersect with annotated genes:"
    wc -l "${ANNOT_INT_LIST}"
    echo "miniprot hits after first filtering:"
    awk '$3=="mRNA"' "${NOCROSS}" | wc -l


    ### check which miniprot IDs overlap each other and only keep the best sequence identity from those
    CROSS_INT_LIST="${SPECIES}_annotation_cross_intersecting.txt"
    awk '$3=="mRNA"' "${NOCROSS}" > mp.gff
    bedtools intersect -wa -wb -s -a mp.gff -b mp.gff | awk -F'\t' '{
        match($9,/ID=[^;]+/); id1=substr($9,RSTART+3,RLENGTH-3);
        match($18,/ID=[^;]+/); id2=substr($18,RSTART+3,RLENGTH-3);
        print id1, id2
    }' > "${CROSS_INT_LIST}"
    # get list of the nonself overlaps that are the best alignent of the overlapping group
    OUTLIST_WEAK_OVERLAP="${SPECIES}_annotation_cross_intersecting_noself.txt"
    
    # check how many IDs there are
    # echo "unique lines in mp.gff"
    # awk '$3=="mRNA"{match($9,/ID=[^;]+/); print substr($9,RSTART+3,RLENGTH-3)}' mp.gff | sort -u | wc -l

    python3 /Users/miltr339/work/PhD_code/PhD_chapter2/src/miniprot_filter_cross_intersection.py "${CROSS_INT_LIST}" "${NOCROSS}" "${OUTLIST_WEAK_OVERLAP}"
    # wc -l ${OUTLIST_WEAK_OVERLAP}
    ANMPLICONS=${SPECIES}_miniprot_no_cross_no_self_hits.gff
    gffread "${NOCROSS}" -F --ids "${OUTLIST_WEAK_OVERLAP}" -o "${ANMPLICONS}"

    echo "miniprot hits only unique non-annotated amplicons:"
    awk '$3=="mRNA"' "${ANMPLICONS}" | wc -l
done 


# test
# for SPECIES in  "A_obtectus" "B_siliquastri" "B_varius" "C_chinensis" "C_maculatus" "D_carinulata" "D_sublineata"  ; do echo $SPECIES; wc -l "${SPECIES}_annotation_non-self_intersecting.txt" ;  ; grep "mRNA" "${SPECIES}_miniprot_no_cross_hits.txt"|wc -l ; echo "---" ;done