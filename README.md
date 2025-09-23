# PhD chapter II
This is all code that was used in chapter2 of my PhD thesis

## synteny

I will start with a synteny analysis of the four bruchid species in this analysis and *Tribolium castaneum* as an outgroup. I will do it the same way as Höök & Näsvall [here](https://doi.org/10.1007/s10577-023-09713-z) with protein synteny via MCScanX and visualization with [SynVisio](https://synvisio.github.io/#/).

### samples

I am using the superscaffolded version of Cmac for now to get a better idea of the big picture. However, since the superscaffolding does not work well for the Y chromosome and just fragments it more, I will use the Kaufmann2023 assembly for the detailed X and Y chromosome related analyses. It is reannotated with BRAKER3 and population-specific RNAseq data (see `bash/reannotate_Kaufmann2023.sh`).

### MCScanX

[MCScanX](https://github.com/wyp1125/MCScanX) is based on protein synteny from blastp searches. The documentation is not super useful beyond the installation instructions, so I was recommended [this workflow](https://www.nature.com/articles/s41596-024-00968-2#Sec29), which I will use for detailed instructions, but for preparing the input files, I am using my own scripts and not their provided ones (esp. the perl scripts have hard-coded aspects that assume ncbi assemblies and annotations, which doesn't apply to my data).

* **make bed files:** I used my script `src/make_bedfile_for_MCScanX.py` which makes all the files according to the format required, including reformatting the contig names, and also creating a lookup table to associate the new contig names with the original ones from the annotations.
*  **blast searches:** I am using blastp according to the above cited workflow on all the species `bash/blastp_for_synteny.sh`  

#### Troubleshooting

* **no alighments generated:** The log file looks like this
    ```
Reading BLAST file and pre-processing
Generating BLAST list
0 matches imported (0 discarded)
0 pairwise comparisons
0 alignments generated
    ```