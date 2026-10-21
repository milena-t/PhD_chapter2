#!/bin/bash -l

#SBATCH -A uppmax2026-1-8
#SBATCH -c 16
#SBATCH -t 50:00:00
#SBATCH -J eggnog_funcannot
#SBATCH -o eggnog_funcannot.log


## Ingo says:
    # Run eggnog-mapper with diamond, running with iterations similar to mmseqs2 with the final iteration on ultra-sensitive settings which should be slightly better than
    # mmseqs2 at s=7.5 (Buchfink, et al. (2025). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature methods)

# Load modules
module load eggnog-mapper/2.1.12-foss-2024a

SPECIES=Cmac


OutDIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/eggnog/${SPECIES}
InFAA=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/protein_data/C_maculatus.faa
GFF=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/native_annotations/C_maculatus.gff

for SPECIES in  A_obtectus B_siliquastri B_varius C_chinensis C_maculatus D_carinulata D_sublineata T_castaneum
do 
    echo $SPECIES
    sbatch  -J eggnog_funcannot_${SPECIES} -o eggnog_funcannot_${SPECIES}.log /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter2/PhD_chapter2/bash/eggnog.sh $SPECIES
done