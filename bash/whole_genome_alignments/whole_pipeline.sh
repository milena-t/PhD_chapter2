#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 16
#SBATCH --mem=200G
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH -J mummer_alignment
#SBATCH -o mummer_alignment.log

# module load bioinfo-tools
module load BEDTools MUMmer/4.0.1-GCCcore-13.3.0
# module load MUMmer/4.0.0rc1

REFERENCE=$1
QUERY=$2

REFERENCE_name="${REFERENCE##*/}"
REFERENCE_name="${REFERENCE_name%.*}"

QUERY_name="${QUERY##*/}"
QUERY_name="${QUERY_name%.*}"

OUT_PREFIX="ref_${REFERENCE_name}_query_${QUERY_name}"

echo "REFERENCE assembly: ${REFERENCE}"
echo "QUERY assembly: ${QUERY}"

echo "outfiles prefix: ${OUT_PREFIX}"

### mummer installed from tarball https://github.com/mummer4/mummer/releases/download/v4.0.1/mummer-4.0.1.tar.gz
# tar -xvzf mummer-4.0.1.tar.gz
# cd mummer-4.0.1
# ./configure --prefix=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/mummer-4.0.1
# make
# make install
MUMMER_PATH=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/mummer-4.0.1

### run mummer alignment
echo " >>>>>>>>>> ALIGNMENT"
echo "nucmer -l 100 --maxmatch -p $OUT_PREFIX -t 16 $REFERENCE $QUERY"
echo "... running nucmer alignment"

# $MUMMER_PATH/nucmer -l 100 --maxmatch -p $OUT_PREFIX -t 16 $REFERENCE $QUERY
NUCMER_OUTFILE="${OUT_PREFIX}.delta"
echo "alignment done! -> ${NUCMER_OUTFILE}"
echo 

echo " >>>>>>>>>> FILTERING"

echo "delta-filter -m -i 50 -l 5000 $NUCMER_OUTFILE > ${NUCMER_OUTFILE}_filtered"
$MUMMER_PATH/delta-filter -m -i 50 -l 2000 $NUCMER_OUTFILE > ${NUCMER_OUTFILE}_filtered
# echo "$MUMMER_PATH/delta-filter -m $NUCMER_OUTFILE > ${NUCMER_OUTFILE}_filtered"
# $MUMMER_PATH/delta-filter -m $NUCMER_OUTFILE > ${NUCMER_OUTFILE}_filtered

FILTERED_ALN=${NUCMER_OUTFILE}_filtered
$MUMMER_PATH/mummerplot -p $OUT_PREFIX -t png $FILTERED_ALN
echo "show-coords -rTH -I 90 ${FILTERED_ALN} > ${FILTERED_ALN}_coords"
$MUMMER_PATH/show-coords -rTH -I 90 $FILTERED_ALN > ${FILTERED_ALN}_coords
COORDS_ALN=${FILTERED_ALN}_coords
echo "filtering done! -> ${COORDS_ALN}"
echo

echo " >>>>>>>>>> BREAKPOINTS"
echo "using the script from the Backström group:"
echo "https://github.com/EBC-butterfly-genomics-team/Leptidea_chromosome_research2022/blob/main/Scripts/Breakpoint_analysis/identify_chr_breakpoins.sh"
SCRIPT_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/bash/whole_genome_alignments
bash $SCRIPT_DIR/identify_chr_breakpoins.sh $COORDS_ALN
echo "breakpoints done!"
