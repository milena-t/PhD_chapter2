## edgeR
## user guide here: https://bioconductor.org/packages//release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# install all dependencies
if (FALSE){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DESeq2") # compiled from source 
  BiocManager::install("edgeR")
  # install tidyverse (ggplot2 is in tidyverse)
  install.packages("tidyverse") 
  install.packages("gtable") 
}
if (TRUE){
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(tibble)
library(stringr)
library(forcats)
library(lubridate)
# library(tidyverse)

library(ggplot2)
library(edgeR)
}

##############
## this is basically identical to chapter3
##############

######################
#### load data
######################
if (TRUE){
  setwd("/Users/miltr339/work/PhD_code/PhD_chapter2/data/DE_analysis")
  Cmac_path <- "/Users/miltr339/work/chapter3/DE_analysis/Cmac_gene_counts_short_headers.txt"
  #Cmac gene lengths
  Cmac_gene_lengths <- read.delim("/Users/miltr339/work/chapter3/DE_analysis/Cmac_gene_lengths.txt", comment.char = "#")

  ## normalized counts outfiles
  Cmac_path_out <- "/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/Cmac_gene_counts_edgeR_normalized.txt"
  ## log-fold-change outfiles
  Cmac_path_out_lfc <- "/Users/miltr339/work/chapter3/DE_analysis/Cmac_gene_counts_short_headers_lfc.txt"
  
  #### metadata relative paths
  Cmac_metadata <- '/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Cmac_SRR_metadata.csv'

  #### do one species of the three
  data_path <- "/Users/miltr339/work/PhD_code/PhD_chapter2/data/DE_analysis"
  out_path_vst <- "/Users/miltr339/work/PhD_code/PhD_chapter2/data/DE_analysis/vst_counts.txt"
  out_path_rpkm <- "/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/Cmac_gene_counts_edgeR_rpkm_length_normalized.txt"
  out_path_lfc <- "/Users/miltr339/work/PhD_code/PhD_chapter2/data/DE_analysis/lfc_counts.txt"
  metadata <- Cmac_metadata
}

######################
##### read data ffrom files and check if metadata matches
######################
if (TRUE){
  ## count data
  count_data <- read.csv(data_path, header = TRUE, sep = "\t",comment.char = "#")
  count_data <- count_data %>% select(-one_of('Chr','Start','End','Strand','Length')) 
  # make the Geneid the row id 
  count_data <- count_data %>% remove_rownames %>% column_to_rownames(var="Geneid")
  ## meta data
  DE_metadata_all <- read.csv(metadata, header = TRUE, sep = ",",comment.char = "#")
  # only include metadata when the samples are in the counts
  meta_data <- DE_metadata_all[which((DE_metadata_all$ID %in% names(count_data))==TRUE), ]
  meta_data <- meta_data %>% remove_rownames %>% column_to_rownames(var="ID")
  # reorder so it's in the same order as the count data columns
  meta_data <- meta_data[names(count_data),]
  # check assumptions
  all(colnames(count_data) %in% rownames(meta_data))
}

#####################
#### read into EdgeR datastructures and do some preprocessing
#####################
if (TRUE){
  # make sample groups that reflect tissue (organ) and sex
  groups <- interaction(meta_data$organ, meta_data$sex, sep="_")
  # read data into data structure
  DE_counts <- DGEList(counts=count_data, group=groups)
  #dim(DE_counts)
  #[1] 14861    12
  
  ## filter out low count genes
  keep <- filterByExpr(DE_counts, group=groups)
  DE_filtered <- DE_counts[keep, , keep.lib.sizes=FALSE]
  DE_filtered$samples$lib.size <- colSums(DE_filtered$counts) # recompute library sizes
  # dim(DE_filtered)
  # [1] 11337    12
  
  ## normalize by library sizes
  DE_normalized <- calcNormFactors(DE_filtered)
}
## export normalized counts to file
if (TRUE){
  log_norm_counts <- cpm(DE_normalized,normalized.lib.sizes = TRUE,log = TRUE,prior.count = 1)
  write.table(log_norm_counts,file = out_path_vst,sep = "\t",quote = FALSE)# ,col.names = NA)
}
## export the gene-length normalized counts file for DC
if (TRUE){
  length_norm_counts <- rpkm(DE_counts, gene.length = Cmac_gene_lengths$Length)
  write.table(length_norm_counts,file = out_path_rpkm,sep = "\t",quote = FALSE)# ,col.names = NA)
}


#####################
#### set up DE analysis
#####################
if (TRUE){
  MM <- model.matrix(~0 + groups) # suppress intercept with '0+'
  rownames(MM) <- colnames(DE_normalized)
  # colnames(MM) <- c("abdomen_female", "abdomen_male", "head_female","head_male") ### check if the order is right before assigning!!!
  colnames(MM) <- c("abdomen_female", "head_female", "abdomen_male","head_male") ### check if the order is right before assigning!!!
  
  # Estimates common, trended and gene-wise dispersion, allowing for a possible trend with average count size 
  DE_dispersion <- estimateDisp(DE_normalized, MM)
  plotBCV(DE_dispersion)
  
  # define contrasts
  my.contrasts <- makeContrasts(
    Sex_a =	abdomen_female - abdomen_male, 
    Sex_h =	head_female - head_male, 
    levels=MM)
  
  # fit linear model
  fit <- glmFit(DE_dispersion, MM)
}


#####################
#### analyze contrasts
#####################
if (TRUE){
  lrt_a <- glmLRT(fit, contrast=my.contrasts[,"Sex_a"])
  summary(de<-decideTestsDGE(lrt_a, p.value=0.05, lfc=1))#, adjust.method="fdr") 
  '1*abdomen_female -1*abdomen_male --> downregulated=male-biased
  Down                               3265
  NotSig                             5279
  Up                                 2793'
  lrt_h <- glmLRT(fit, contrast=my.contrasts[,"Sex_h"])
  summary(de<-decideTestsDGE(lrt_h, p.value=0.05, lfc=1))#, adjust.method="fdr") 
  '1*head_female -1*head_male
  Down                         1331
  NotSig                       9201
  Up                            805'
  
  ### try to reproduce above numbers from topTags()
  # abdomen (correct!)
  if (FALSE){
    SexbAbd <- topTags(lrt_a, n=Inf)
    SexbAbd_dat <-SexbAbd@.Data[[1]] 
    length(SexbAbd_dat$logFC) # 11337 -> full correct dataset
    SexbAbd_sig <- SexbAbd_dat[SexbAbd_dat$FDR < 0.05, ]
    length(SexbAbd_sig$logFC) # 9353
    
    
    SexbAbd <- topTags(lrt_a, n=Inf, p.value=0.05)	
    SexbAbd_dat <-SexbAbd@.Data[[1]] 
    length(SexbAbd_dat$logFC) # 9353
    SexbAbd_up <- SexbAbd_dat[SexbAbd_dat$logFC > 1, ] # 2793
    SexbAbd_do <- SexbAbd_dat[SexbAbd_dat$logFC < -1, ] # 3265
  }
  
  ## save LFC
  ## LFC is female-male !!
  SexbAbd <- topTags(lrt_a, n=Inf)	
  write.table(SexbAbd,file = "/Users/miltr339/work/PhD_code/PhD_chapter2/data/DE_analysis/Cmac_sex_DE_abdomen_edgeR.txt",sep = "\t",quote = FALSE)# ,col.names = NA)
  SexbHT <- topTags(lrt_h, n=Inf)
  write.table(SexbHT,file = "/Users/miltr339/work/PhD_code/PhD_chapter2/data/DE_analysis/Cmac_sex_DE_head_thorax_edgeR.txt",sep = "\t",quote = FALSE)# ,col.names = NA)
}



