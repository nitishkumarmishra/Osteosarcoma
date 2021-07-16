setwd("D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma")
load("D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/TARGET-OS-RNAseq.RData")


countdata <- OS_HTSeq_Counts
countdata <- log2(countdata + 1)
cnts <- countdata[rowSums(countdata==0)< ncol(countdata)*0.9,] ## Remove all gene which have 25% zero's
keep <- rowSums(edgeR::cpm(cnts)>1) >= ncol(cnts)*0.50 #### Select only genes which have have CPM > 1 for >=50% samples
cnts <- cnts[keep,]; rm(keep)
cnts <- as.data.frame(cnts)

gencode.v34.annotation <- readRDS("gencode.v34.annotation.rds")
gencode.v34.annotation$ENSG <- gsub("\\.[1-9]", "", gencode.v34.annotation$ENSG)

##############################################################################
##############################################################################
## This part is to generate gencode.v34.gtf. I already save an loaded by using previous line
# gtf <- rtracklayer::import('gencode.v34.annotation.gtf.gz')
# gencode.v34.gtf <- as.data.frame(gtf)
# 
# library(dplyr)
# library(tidyverse)
# 
# gencode.v34.gtf.selected <- gencode.v34.gtf %>%
#   filter(type=="gene") %>%
#   rename(Chr=seqnames, Type=type,ENSG=gene_id, GeneType=gene_type, Symbol=gene_name) %>%
#   select(Chr, Type, GeneType, ENSG, Symbol) 
# 
# saveRDS(gencode.v34.gtf.selected, file = "gencode.v34.annotation.rds")
# ##############################################################################
##############################################################################
library(dplyr)
library(tidyverse)
cnts$ENSG <- rownames(cnts)
cnts <- cnts %>%
  inner_join(gencode.v34.annotation, by = "ENSG")
