# PhD chapter II
This is all code that was used in chapter2 of my PhD thesis

# synteny

I will start with a synteny analysis of the four bruchid species in this analysis and *Tribolium castaneum* as an outgroup. I will do it the same way as Höök & Näsvall [here](https://doi.org/10.1007/s10577-023-09713-z) with protein synteny via MCScanX and visualization with [SynVisio](https://synvisio.github.io/#/).

## samples

I am using the superscaffolded version of Cmac for now to get a better idea of the big picture. However, since the superscaffolding does not work well for the Y chromosome and just fragments it more, I will use the Kaufmann2023 assembly and annotation for the detailed X and Y chromosome related analyses. The annotation is not ideal in some aspect, since it estimates a very high number of genes and has a high proportion of single exon genes which I hypothesize are caused by some default values in the BRAKER2 and TSEBRA pipeline being not ideal for these kinds of large beetle genomes. I have attempted to reannotate the assembly with BRAKER3 and the same population-specific RNAseq data as the original annotation (see `bash/reannotate_Kaufmann2023.sh`), but this struggles to identify the TOR copy number variation on the Y chromosome (in a way that it does not when RNAseq data is excluded), therefore I have decided to continue using the Kaufmann 2023 version of the annotation. See detailed report on me figuring it out in [chapter 4](https://github.com/milena-t/PhD_chapter4/blob/main/mTOR_annotation/mTOR_notes.md), especially the figure at the end.

## MCScanX

[MCScanX](https://github.com/wyp1125/MCScanX) is based on protein synteny from blastp searches. The documentation is not super useful beyond the installation instructions, so I was recommended [this workflow](https://www.nature.com/articles/s41596-024-00968-2#Sec29), which I will use for detailed instructions, but for preparing the input files, I am using my own scripts and not their provided ones (esp. the perl scripts have hard-coded aspects that assume ncbi assemblies and annotations, which doesn't apply to my data).

* **make bed files:** I used my script `src/make_bedfile_for_MCScanX.py` which makes all the files according to the format required, including reformatting the contig names, and also creating a lookup table to associate the new contig names with the original ones from the annotations.
*  **blast searches:** I am using blastp according to the above cited workflow on all the species `bash/blastp_for_synteny.sh`  

### Troubleshooting

* **no alighments generated:** The log file reads no blast alignments. For me, this was caused by the file extension of the annotation file. I named it `.bed` but it should be `.gff` (because it is somehow hardcoded?)

## Results plotted with SynVisio

I have visualized the results with [SynVisio](https://synvisio.github.io/#/).
The order in the species in the below plot from top to bottom is:

* *Acanthoscelides obtectus* (ao)
* *Callosobruchus maculatus* (cm)
* *Bruchidius siliquastri* (bs)
* *Tribolium castaneum* (tc)
  
*C. chinensis* was included in the analysis but not the plot because the assembly is so fragmented it is not helpful to visualize here.

![synteny plot](data/images/synvisio_plot_Cmac_on_top.png)

The sex chromosomes are these (chromosomes with no syntenic regions are excluded from the plot)

* ***B. siliquastri***: X: `bs9`, Y: `bs8`
  
TODO others

* ***A. obtectus***:
* ***C. maculatus***:
* ***T. castaneum***:

see the [alternative synteny plot](data/images/synvisio_plot.png) for a different order of the bruchid species. Their tree topology is `((Cmac,Bsil),Aobt)`, but there is a bunch of rearrangements between all of them, and *C. maculatus* is not visibly more similar to *B. siliquastri* than to *A. obtectus*.

# Translocation analysis

I aim to investigate the processes by which the genes migrate to the sex chromosomes resulting in the patterns we observe today. Orthologs, including gametologs, are identified with orthofinder between all species, and with the help of the known phylogeny I will trace back when and how genes on the Y chromosome originated.

## Notes

* [Martinez-Pacheco 2020](https://academic.oup.com/gbe/article/12/11/2015/5892261) 
  * double check y linked gene presence by `blastn` against the assembly.
  * *"the best BlastN match (usually around 92–95% identity over the entire sequence) onto the annotated X chromosome of the reference genomes was considered the X gametologs"*
  * If the X gametolog is missing in the annotation, the sequene from the transcriptome was used instead
  * [Marques 2005](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0030357) for identifying retrogenes. blast proteins against assembly, merge nearby matches, *"query and target sequences had >50% similarity on the amino acid level and over >80% of their length* \[are\] *shared"*, verify absence of introns. with some `paml` stuff, they identified the ancestral gametolog that all retrogenes originate from (useful for FastX?)
* [Peneder 2017](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.3278)
  * 

## orthofinder results
