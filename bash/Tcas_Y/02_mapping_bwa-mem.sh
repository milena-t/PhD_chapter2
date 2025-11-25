#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 20:00:00
#SBATCH -J mapping_bwa
#SBATCH -o mapping_bwa.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se


module load bioinfo-tools  python3  bwa/0.7.17 samtools/1.17 bwa-mem2/2.2.1-20211213-edc703f

ref_dir=/proj/naiss2023-6-65/Milena/chapter2/Tribolium_poolseq

TCAS_index=T_castaneum_index_for_bwa

### index reference before mapping 
bwa index -p $TCAS_index "${ref_dir}/T_castaneum_assembly.fna"

# Directory with reads
reads_dir="/proj/naiss2023-6-65/Milena/chapter2/Tribolium_poolseq/raw_data/trimmed_fastp"
mapped_dir="/proj/naiss2023-6-65/Milena/chapter2/Tribolium_poolseq/mapped_data"

for SRR_num in  SRR23732240 SRR23732241 SRR23732242 SRR23732243 SRR23732244 SRR23732245
do
    sbatch -J "${SRR_num}_mapping" -o "${SRR_num}_mapping.out" /proj/naiss2023-6-65/Milena/chapter2/PhD_chapter2/bash/Tcas_Y/021_mapping_bwa_mem_ind.sh $SRR_num
done