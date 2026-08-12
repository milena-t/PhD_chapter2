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
* I adjust the braker4 flags as follows 
  * `--keep-going` even if one sample fails continue with the other ones
  * `--snakefile /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/software_install/BRAKER4/Snakefile`
  * `--executor slurm`
  * `--default-resources slurm_partition=pelle mem_mb=120000` pelle partition
  * `--cores 20`
  * `--jobs 50` max. jobs at one time
  * `--use-singularity`
  * `--singularity-prefix .singularity_cache`
  * `--singularity-args "-B /home -B /scratch"`
  * `--latency-wait 120`
  * `--restart-times 3`

run on an interactive node for the first time when it makes the singularity image, and then cancel. login node for all subsequent runs, as is recommended by snakemake.

