The new biopython update breaks the braker3 singularity container, so I will try to use braker4.
This tutorial: https://github.com/Gaius-Augustus/BRAKER4/blob/main/MIGRATING_FROM_BRAKER3.md

* the git repository is here `/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/software_install/BRAKER4`
* I have these command line flags in BRAKER3 
  * `--genome=${ASSEMBLY_MASKED}`
  * `--prot_seq $PROTEIN_DATA`
  * `--rnaseq_sets_ids=$SRA_IDS`
  * `--threads 20`
  * `--GENEMARK_PATH=${ETP}/gmes`
  * `--AUGUSTUS_CONFIG_PATH=${wd}/AUGUSTUS_config`
  * `--useexisting`

* I adjust `config.ini`: 20-core SLURM node with 120 GB RAM and a 72-hour wall time.