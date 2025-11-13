#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
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

SRR_num=$1
r1="${output_dir}/${SRR_num}_1_trimmed.fastq.gz"
r2="${output_dir}/${SRR_num}_2_trimmed.fastq.gz"
echo "======================>> Running bwa-mem on $sample ..."
bwa mem -t 20 -P $TCAS_index $r1 $r2 | samtools view -u | samtools sort -o "${mapped_dir}/${SRR_num}.bam"
