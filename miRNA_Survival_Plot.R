setwd("D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/TARGET-OS/")
load("TARGET_OS_Survival.RData")
rm(list=setdiff(ls(), c("OS.Clinical", "OS.Clinical.TARGET", "dataset.protein.coding")))

#######
##miRNA expression of TARGET-OS
## With normalization
genedat = "ftp://caftpd.nci.nih.gov/pub/OCG-DCC/TARGET/OS/miRNA_pcr/L2/TARGET_OS_miRNA_level2_dCt.txt"
OS.miRNA = readr::read_tsv(genedat)

library(tidyverse)
#OS.miRNA <- OS.miRNA
OS.miRNA <- OS.miRNA[!duplicated(OS.miRNA$miRNA),]
OS.miRNA <- OS.miRNA %>%
  tibble::column_to_rownames("miRNA")

## sapply(strsplit(rownames(miRNA), "-"), function(x) paste(x[1:3], collapse = "-")) ## Take first three 
## sapply(strsplit(rownames(miRNA), "-"), function(x) paste(head(x, n=-1), collapse = "-")) ## Take all except last


source("plot_miRNA_Surv_RNAseq_function.R")
#unlink("RNAseq_Surv/*") 
#results <- plot_surv.Protein(dir = "RNAseq_Surv", clinical_patient = OS.Clinical.TARGET, dataGE = dataset.protein.coding, Genelist = rownames(dataset.protein.coding), Survresult = FALSE, Median = TRUE, p.cut = 0.05)

unlink("RNAseq_miRNA_Surv_P_0.01/*") 
results_miRNA.P.0.01 <- plot_surv.Protein(dir = "RNAseq_miRNA_Surv_P_0.01", clinical_patient = OS.Clinical.TARGET, dataGE = OS.miRNA, Genelist = rownames(OS.miRNA), Survresult = TRUE, Median = TRUE, p.cut = 0.01)

unlink("RNAseq_miRNA_Surv_P_0.05/*") 
results_miRNA.P.0.05 <- plot_surv.Protein(dir = "RNAseq_miRNA_Surv_P_0.05", clinical_patient = OS.Clinical.TARGET, dataGE = OS.miRNA, Genelist = rownames(OS.miRNA), Survresult = TRUE, Median = TRUE, p.cut = 0.05)

save.image("miRNA_Survival.RData")
