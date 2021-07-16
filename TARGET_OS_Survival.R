## ---------------------------
library(TCGAbiolinks)
library("tidyr"); library(tidyverse); library(dplyr)
library("tidylog", warn.conflicts = FALSE)
library(matrixStats)
######################################################
setwd("D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/TARGET-OS/")

########## Download clinical data ####################
OS.Clinical <- GDCquery_clinic("TARGET-OS")

########### HTSeq - Counts download ##################
### We need TARGET-OS for expression and methylation data, but only OS for MAF files
query<- GDCquery(project = "TARGET-OS", data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "HTSeq - Counts")
GDCdownload(query)
OS.HTSeq.Counts <- GDCprepare(query = query, summarizedExperiment = FALSE)
OS.HTSeq.Counts <- tibble::column_to_rownames(OS.HTSeq.Counts, var = "X1")

########## HTSeq - FPKM download #####################
query<- GDCquery(project = "TARGET-OS", data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "HTSeq - FPKM")
GDCdownload(query)
OS.HTSeq.FPKM <- GDCprepare(query = query, summarizedExperiment = FALSE)
OS.HTSeq.FPKM <- tibble::column_to_rownames(OS.HTSeq.FPKM, var = "X1")

######### HTSeq - FPKM-UQ download ###################
query<- GDCquery(project = "TARGET-OS", data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query)
OS.FPKM.UQ <- GDCprepare(query = query, summarizedExperiment = FALSE)
OS.FPKM.UQ <- tibble::column_to_rownames(OS.FPKM.UQ, var = "X1")

## Read lastest clinical file downloaded from the TARGET
## https://target-data.nci.nih.gov/Public/OS/clinical/harmonized/
OS.Clinical.Discovery.2018 <- readxl::read_excel("../Literature/TARGET_OS_ClinicalData_Discovery_20181009.xlsx")
OS.Clinical.Validation.2019 <- readxl::read_excel("../Literature/TARGET_OS_ClinicalData_Validation_20190805.xlsx")
OS.Clinical.TARGET <- rbind(OS.Clinical.Discovery.2018, OS.Clinical.Validation.2019)

##############################################################################
## This part is to generate gencode.v22.gtf. I already save an loaded by using previous line
gtf <- rtracklayer::import('../GDC Osteosarcoma/gencode.v22.annotation.gtf.gz')
gencode.v22.gtf <- as.data.frame(gtf)

gencode.v22.gtf.TCGA.selected <- gencode.v22.gtf %>%
  filter(type=="gene") %>%
  rename(Chr=seqnames, Type=type,ENSG=gene_id, GeneType=gene_type, Symbol=gene_name, Status=gene_status) %>%
  select(Chr, Type, GeneType, ENSG, Symbol)

count <- OS.HTSeq.FPKM
count <- count[rowSums(count==0)< ncol(count)*0.20,] ## Remove all gene which have 25% zero's
keep <- rowSums(edgeR::cpm(count)>1) >= ncol(count)*0.10 #### Select only genes which have have CPM > 1 for >=50% samples
count <- count[keep,]
count <- log2(OS.HTSeq.Counts[rownames(count),]+1)

EnsembleID <- rownames(count)
count <- count %>%
  mutate(ENSG=EnsembleID)

dataset.protein.coding <- count %>%
  inner_join(gencode.v22.gtf.TCGA.selected, by = "ENSG") %>%
  filter(GeneType %in% c("protein_coding"))
dataset.protein.coding <- dataset.protein.coding[!duplicated(dataset.protein.coding$Symbol),] ##Remove duplicate genes
rownames(dataset.protein.coding) <- dataset.protein.coding$Symbol

dataset.lincRNA <- count %>%
  inner_join(gencode.v22.gtf.TCGA.selected, by = "ENSG") %>%
  filter(GeneType %in% c("lincRNA"))
dataset.lincRNA <- dataset.lincRNA[!duplicated(dataset.lincRNA$Symbol),]#Remove duplicate lincRNAs
rownames(dataset.lincRNA) <- dataset.lincRNA$Symbol

OS.Clinical.TARGET <- as.data.frame(OS.Clinical.TARGET)
OS.Clinical.TARGET <- OS.Clinical.TARGET[!duplicated(OS.Clinical.TARGET$`TARGET USI`), ]

##########################################
save.image("TARGET_OS_Survival.RData")
