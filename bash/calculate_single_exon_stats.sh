#!/bin/sh
#SBATCH -A naiss2024-5-135
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 20:00
#SBATCH -J single_exon_stats
#SBATCH -o single_exon_stats.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# module load bioinfo-tools python3

# all orthoDB annotations

# for SPECIES in  A_obtectus B_siliquastri B_varius C_chinensis C_maculatus D_carinulata D_sublineata
# do
#     python3 /Users/miltr339/work/PhD_code/PhD_chapter2/src/calculate_single_exon_stats.py \
#     "/Users/miltr339/work/chapter2/native_annotations/${SPECIES}.gff" \
#     True > "/Users/miltr339/work/PhD_code/PhD_chapter2/data/single_exon_stats/${SPECIES}_single_exon_stats_with_transcript_list.txt"
#     echo "-----> done native annot ${SPECIES}"
# done

for SPECIES in  A_obtectus B_siliquastri B_varius C_chinensis C_maculatus D_carinulata D_sublineata
do
    python3 /Users/miltr339/work/PhD_code/PhD_chapter2/src/calculate_single_exon_stats.py \
    "/Users/miltr339/work/chapter2/miniprot_annot/${SPECIES}_miniprot_no_cross_hits.gff" \
    True True > "/Users/miltr339/work/PhD_code/PhD_chapter2/data/single_exon_stats/${SPECIES}_miniprot_single_exon_stats_with_transcript_list_no_cross_hits.txt"
    echo "-----> done miniprot ${SPECIES}"
done