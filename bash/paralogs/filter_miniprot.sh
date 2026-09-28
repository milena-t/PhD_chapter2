#!/bin/bash -l

MINI_DIR=/Users/miltr339/work/chapter2/miniprot_annot
ANN_DIR=/Users/miltr339/work/chapter2/native_annotations

## filter for only mapped regions that aren't annotated with non-query genes

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

    ## intersect by position, keep IDs where a mapped miniprot region intersects with annotated genes except non-self hits
    ANNOT_INT_LIST="${SPECIES}_annotation_self_intersecting.txt"
    bedtools intersect -a mp_mRNA.gff -b genes.gff -wa -wb \
    | awk -F'\t' '{
        match($9,/(^|;)ID=[^;]+/);      id=substr($9,RSTART,RLENGTH);  sub(/^;?ID=/,"",id);
        match($9,/Target=[^ ;]+/);      t=substr($9,RSTART+7,RLENGTH-7);
        match($18,/(^|;)ID=[^;]+/);     g=substr($18,RSTART,RLENGTH);  sub(/^;?ID=/,"",g);
        if (t!=g) print id
    }' | sort -u > "${ANNOT_INT_LIST}"

    ## remove all IDs in the list generated above
    # this does not include the #PAF lines because gffread does not like them
    NOCROSS=${SPECIES}_miniprot_no_cross_hits.gff
    grep -v '^##PAF' "${MINIPROT}" | gffread - -F --nids "${ANNOT_INT_LIST}" -o "${NOCROSS}"

    echo "   alignent IDs that intersect with other (non-query) annotated genes:"
    wc -l "${ANNOT_INT_LIST}"
    echo "miniprot hits after first filtering:"
    awk '$3=="mRNA"' "${NOCROSS}" | wc -l


    ### check which miniprot IDs overlap each other and only keep the best sequence identity from those
    ## TODO this does run really long
    MINI_SELF_LIST="${SPECIES}_miniprot_self_intersecting.txt"
    awk -F'\t' -v OFS='\t' '$3=="mRNA"{
        match($9,/ID=[^;]+/);       id=substr($9,RSTART+3,RLENGTH-3);
        match($9,/Identity=[^;]+/); idt=substr($9,RSTART+9,RLENGTH-9);
        print id,$1,$4,$5,$7,idt,$6}' "${NOCROSS}" > "${MINI_SELF_LIST}"
    NOSELF=${SPECIES}_miniprot_no_cross_no_self_hits.gff
    echo "   alignent IDs that intersect with other alignment:"
    wc -l "${ANNOT_INT_LIST}"
    
    gffread "${NOCROSS}" -F --nids "${MINI_SELF_LIST}" -o "${NOSELF}"
    echo "miniprot hits after second filtering:"
    if [ "${FILTSTR}" = "transcript" ] ; then
        awk '$3=="transcript"' "${NOESELF}" | wc -l
    else 
        awk '$3=="mRNA"' "${NOESELF}" | wc -l
    fi
done 


# test
# for SPECIES in  "A_obtectus" "B_siliquastri" "B_varius" "C_chinensis" "C_maculatus" "D_carinulata" "D_sublineata"  ; do echo $SPECIES; wc -l "${SPECIES}_annotation_non-self_intersecting.txt" ;  ; grep "mRNA" "${SPECIES}_miniprot_no_cross_hits.txt"|wc -l ; echo "---" ;done