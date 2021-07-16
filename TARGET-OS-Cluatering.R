## ---------------------------
## Script name: TARGET_OS_Analysis.R
## Purpose of script: TARGET-OS data download and Consensus clustering analysis
## Author: Dr. Nitish Kumar Mishra
## Date Created: 2020-07-04
## Copyright (c) Nitish Kumar Mishra, 2020
## Email: nitish.mishra@unmc.edu
## ---------------------------
library(TCGAbiolinks)
library("tidyr"); library(tidyverse); library(dplyr)
library("tidylog", warn.conflicts = FALSE)
library(matrixStats); library(maftools)
library(SNFtool); library(TargetOsteoAnalysis)
library(iClusterPlus); library(GenomicRanges)
library(gplots); library(lattice); library(cluster)
library(ConsensusClusterPlus)
######################################################
dir.create(file.path("TARGET-OS"), showWarnings = FALSE)
setwd("D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/TARGET-OS/")

########## Download clinical data ####################
OS.Clinical <- GDCquery_clinic("TARGET-OS")

########### HTSeq - Counts download ##################
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
## This part is to load gencode.v22.gtf
gencode.v22.gtf.TCGA.selected <- readRDS(file = "../GDC Osteosarcoma/gencode.v22.annotation.TCGA.rds")
Osteo <- read.maf(maf = "../GDC Osteosarcoma/Combined_All_MuTect2.aliquot.maf", clinicalData = NULL)

#########################################################
################### TARGET-OS CNV #######################
variants = variant_calls()
#OS.GISTIC <- target_gistic_se()
gistic_file = system.file(package='TargetOsteoAnalysis','extdata/all_data_by_genes.txt.bz2')
OS.GISTIC = read_tsv(gistic_file)
OS.Purity = target_purity()

fname = system.file("extdata/summary.pureCN.cns.calls.tsv.gz", package = "TargetOsteoAnalysis")
OS.CNVkit.Purity = read_tsv(fname, col_names = TRUE)
############## miRNA expression of TARGET-OS ############
#genedat = "ftp://caftpd.nci.nih.gov/pub/OCG-DCC/TARGET/OS/miRNA_pcr/L2/TARGET_OS_miRNA_level2_dCt.txt"
genedat = "https://target-data.nci.nih.gov/Public/OS/miRNA_pcr/L2/TARGET_OS_miRNA_level2_dCt.txt"
OS.miRNA = readr::read_tsv(genedat)
## without normalization miRNA 
genedat = "https://target-data.nci.nih.gov/Public/OS/miRNA_pcr/L2/TARGET_OS_miRNA_level2_Ct.txt"
OS.miRNA1 = readr::read_tsv(genedat)
#########################################################
################### DNA Methylation #####################
OS.Meth <- readRDS("../OsteoMeth.rds")
hm450.anno <-read.table("../hm450.hg38.manifest.tsv.gz",header=T,sep="\t",na.strings="",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)
rownames(hm450.anno) <- hm450.anno$probeID
hm450.anno.unmasked <- hm450.anno %>% 
  filter(MASK_general %in% "FALSE") %>%
  filter(MASK_snp5_common %in% "FALSE")
masked <- which(hm450.anno.unmasked$chrm_A %in% c("chrX", "chrY", "chrM"))
hm450.anno.unmasked <- hm450.anno.unmasked[-masked,]

Meth <- OS.Meth[rownames(OS.Meth)%in%rownames(hm450.anno.unmasked),]
Meth.OS <- Meth[,grep("TARGET-", colnames(Meth))]
Meth.OS <- t(Meth.OS)
###########################################################
################# iClusterPlus ############################
Segment <- readr::read_tsv("../GenePattern GISTIC/TARGET_OS_L3_Segmentation.txt ", col_names = TRUE) %>%
  mutate(Chromosome = gsub("chrX", "chr23", Chromosome)) %>%
  mutate(Chromosome = gsub("chrY", "chr24", Chromosome)) %>%
  filter(Chromosome != "chrMT") %>% 
  mutate(Chromosome = gsub("chr", "", Chromosome)) %>%
  mutate(Chromosome = as.numeric(Chromosome))

TCGA.CNV.bypass <- read.csv("../GenePattern GISTIC/CNV.hg19.bypos.111213.txt", header = TRUE, sep = "\t")
TCGA.CNV.bypass$Chromosome <- paste0("chr", TCGA.CNV.bypass$Chromosome)

seg <- as.data.frame(Segment)
OS.cn=CNregions(seg=seg,epsilon=0.02,adaptive=FALSE,rmCNV=TRUE,
                cnv=TCGA.CNV.bypass[,2:4],
                frac.overlap=0.5, rmSmallseg=TRUE,nProbes=5) # Total 6019 regions

Segment <- Segment %>%
  mutate(Source = substr(Segment$SampleName, 11, 16))

Segment.Purity <- left_join(Segment, OS.Purity, by=c("Source")) %>%
  select(!c(Source, dataset)) %>%
  mutate("Segment Purity" = if_else(!is.na(Purity),((`Segment Mean`/Purity) - 2*(1-Purity)/(Purity*Ploidy)),`Segment Mean`))

Segment.Purity.selected <- Segment.Purity %>%
  select(c(SampleName, Chromosome, Start, End, `Num Probes`, `Segment Purity`)) %>%
  dplyr::rename("Segment Mean"="Segment Purity")

seg.purity <- as.data.frame(Segment.Purity.selected)
OS.cn.purity <- CNregions(seg=seg.purity,epsilon=0.05,adaptive=TRUE,rmCNV=TRUE,
                          cnv=TCGA.CNV.bypass[,2:4],
                          frac.overlap=0.5, rmSmallseg=TRUE,nProbes=5) # Total 2122 regions

OS.cn.purity.1 <- CNregions(seg=seg.purity,epsilon=0.04,adaptive=FALSE,rmCNV=TRUE,
                            cnv=TCGA.CNV.bypass[,2:4],
                            frac.overlap=0.6, rmSmallseg=TRUE,nProbes=5) # Total 5996 regions

####################################################################
Osteo.MuTect.GISTIC <- read.maf(maf = "../GDC Osteosarcoma/Combined_All_MuTect2.aliquot.maf", clinicalData = NULL,
                                gisticAllLesionsFile="../GenePattern GISTIC/269276/all_lesions.conf_90.txt",
                                gisticAmpGenesFile= "../GenePattern GISTIC/269276/amp_genes.conf_90.txt",
                                gisticDelGenesFile = "../GenePattern GISTIC/269276/del_genes.conf_90.txt",
                                gisticScoresFile = "../GenePattern GISTIC/269276/scores.gistic")

MuTect.matrix <- mutCountMatrix(maf = Osteo.MuTect.GISTIC)
## Remove genes which has too few mutations for clustering (less than 5%)
MuTect.matrix <- MuTect.matrix[apply(MuTect.matrix==0, 1, sum) <= ncol(MuTect.matrix)*0.95,]
MuTect.matrix.Binary <- ifelse(MuTect.matrix>0, 1, 0)
MuTect.matrix.Binary <- t(MuTect.matrix.Binary)
mut.rate <- apply(MuTect.matrix.Binary,2,mean)
MuTect.matrix.Binary2 = MuTect.matrix.Binary[,which(mut.rate>0.2)]

###################################################################
count <- OS.HTSeq.FPKM
count <- count[rowSums(count==0)< ncol(count)*0.20,] ## Remove all gene which have 25% zero's
keep <- rowSums(edgeR::cpm(count)>1) >= ncol(count)*0.10 #### Select only genes which have have CPM > 1 for >=50% samples
count <- count[keep,]
#count <- log2(count+1)
############### gene median center normalization
############### https://www.sciencedirect.com/science/article/pii/S1535610817300016?via%3Dihub
#https://www.cell.com/cms/10.1016/j.ccell.2017.01.001/attachment/d7e976b2-aac3-4140-a993-5d212979a499/mmc1.pdf
count <- count-apply(count, 1, median) 
EnsembleID <- rownames(count)
count <- count %>%
  mutate(ENSG=EnsembleID)
dataset.protein.coding <- count %>%
  inner_join(gencode.v22.gtf.TCGA.selected, by = "ENSG") %>%
  filter(GeneType %in% c("protein_coding"))
rownames(dataset.protein.coding) <- dataset.protein.coding$ENSG
dataset.protein.coding[,c("ENSG", "Chr", "Type", "GeneType", "Symbol")] <- NULL
#The Similar Network Fusion (SNF) method was run on primary tumor samples using both gene expression and DNA methylation
#data (Wang et al., 2014). The SNF method does not require any prior feature selection so we used the full matrix of gene expression
#(21,641 genes) and the full matrix of methylation data (logitB values, 321,174 probes).

miRNA <- OS.miRNA
miRNA <- miRNA[!duplicated(miRNA$miRNA),]
miRNA <- miRNA %>%
  tibble::column_to_rownames("miRNA")
miRNA <- t(miRNA)
expression <- t(dataset.protein.coding)
rownames(expression) <- substr(rownames(expression), 1, 16)
mutation <- MuTect.matrix.Binary2
rownames(mutation) <- substr(rownames(mutation), 1, 16)
copyNumber <- OS.cn.purity.1
rownames(copyNumber) <- substr(rownames(copyNumber), 1, 16)
methylation <- Meth.OS
rownames(methylation) <- substr(rownames(methylation), 1, 16)

GISTIC <- read.csv("../GenePattern GISTIC/269276/all_lesions.conf_90.txt", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
rownames(GISTIC) <- GISTIC$`Unique Name`
GISTIC <- GISTIC[,grep("TARGET-", colnames(GISTIC))]
GISTIC <- t(GISTIC)
rownames(GISTIC) <- substr(rownames(GISTIC), 1, 16)

#list_of_data = list(methylation , expression, miRNA, mutation, copyNumber, GISTIC)
list_of_data = list(methylation , expression, miRNA, mutation, GISTIC)
common_names = Reduce(intersect, lapply(list_of_data, row.names))

methylation1 <- methylation[common_names,]
expression1 <- expression[common_names,]
###################################################################
############################## SNFTools ###########################
## First, set all the parameters:
K = 20;		# number of neighbors, usually (10~30)
alpha = 0.6;  	# hyperparameter, usually (0.3~0.8)
T = 50; 	# Number of Iterations, usually (10~20)
## In data rows are samples and columns are gene/probes 
Data1 = lumi::beta2m(methylation1)
Data2 = expression1
#Data1 = standardNormalization(methylation1);
#Data2 = standardNormalization(expression1);

## Calculate the pair-wise distance; 
## If the data is continuous, we recommend to use the function "dist2" as follows 
Dist1 = (dist2(as.matrix(Data1),as.matrix(Data1)))^(1/2)
Dist2 = (dist2(as.matrix(Data2),as.matrix(Data2)))^(1/2)

## next, construct similarity graphs
W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)

## These similarity graphs have complementary information about clusters.
#displayClusters(W1, truelabel);
#displayClusters(W2, truelabel);

## next, we fuse all the graphs, then the overall matrix can be computed by similarity network fusion(SNF):
W = SNF(list(W1,W2), K, T)
## With this unified graph W of size n x n, I can do either spectral clustering or Kernel NMF. 
C = 5 								# number of clusters
group = spectralClustering(W,C); 	# the final subtypes information
displayClusters(W, group)
labels = spectralClustering(W, C)
displayClustersWithHeatmap(W, group = group)
## You can get cluster labels for each data point by spectral clustering
plot(Data1, col=labels, main='Data type 1')
plot(Data2, col=labels, main='Data type 2')
estimationResult = estimateNumberOfClustersGivenGraph(W, 2:5);
###################################################################
######################## clustering on data #######################
count <- OS.HTSeq.FPKM
colnames(count) <- substr(colnames(count), 1, 16)
count <- count[,common_names]
count <- count[rowSums(count==0)< ncol(count)*0.20,] ## Remove all gene which have 25% zero's
keep <- rowSums(edgeR::cpm(count)>1) >= ncol(count)*0.10 #### Select only genes which have have CPM > 1 for >=50% samples
count <- count[keep,]
count <- log2(count+1)
#count = sweep(count,1, apply(count,1,median,na.rm=T)) #ConsensusClusterPlus Reference Manual
count <- count-apply(count, 1, median) # sweep and this give same value

EnsembleID <- rownames(count)
count <- count %>%
  #rownames_to_column("ENSG") %>%
  mutate(Mean= rowMeans(.), Median=rowMedians(as.matrix(.)) ,Var=rowVars(as.matrix(.)), 
         stdev=rowSds(as.matrix(.)), MAD=rowMads(as.matrix(.), method = "median")) %>%
  mutate(ENSG=EnsembleID)

Num_Gene=3000
dataset.protein.coding <- count %>%
  inner_join(gencode.v22.gtf.TCGA.selected, by = "ENSG") %>%
  filter(GeneType %in% c("protein_coding")) %>%
  slice_max(MAD, n=Num_Gene)

Num_Gene=500
dataset.lincRNA <- count %>%
  inner_join(gencode.v22.gtf.TCGA.selected, by = "ENSG") %>%
  filter(GeneType %in% c("lincRNA")) %>%
  slice_max(MAD, n=Num_Gene)

############### concensusClusterPlus ####################
## On Median Center Normalizes FPKM; MAD based, https://bioconductor.org/packages/release/bioc/vignettes/ConsensusClusterPlus/inst/doc/ConsensusClusterPlus.pdf
d.FPKM.coding <- as.matrix(dataset.protein.coding[,grep("TARGET-", colnames(dataset.protein.coding))])

dir.create("12345"); title="12345"
results.coding = ConsensusClusterPlus(d.FPKM.coding,maxK=10,reps=1000,pItem=1,pFeature=0.9, title=title,
                                      clusterAlg="hc",distance="pearson", seed=12345,plot="png")

d.FPKM.lincRNA <- as.matrix(dataset.lincRNA[,1:88])
dir.create("12346"); title="12346"
results.lincRNA = ConsensusClusterPlus(d.FPKM.lincRNA,maxK=10,reps=1000,pItem=1,pFeature=0.9, title="lincRNA",
                                       clusterAlg="hc",distance="pearson", seed=12345,plot="png")


dir.create("12347"); title="12347"
results.coding = ConsensusClusterPlus(t(OS.cn.purity.1),maxK=10,reps=1000,pItem=1,pFeature=0.9, title="CNV",
                                      clusterAlg="hc",distance="pearson", seed=12345,plot="png")


###################################################################
######################### CNV clustering ##########################
#setwd("D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/TARGET-OS/")
CNV_by_gene <- readr::read_tsv("../GenePattern GISTIC/GISTIC Purity Ploidy Q 0.25/OS_Purity_Ploidy_q_default.all_data_by_genes.txt") %>%
  column_to_rownames("Gene Symbol") %>%
  select_if(grepl("TARGET", names(.))) %>%
  distinct()

CNV_by_threshold <- readr::read_tsv("../GenePattern GISTIC/GISTIC Purity Ploidy Q 0.25/OS_Purity_Ploidy_q_default.all_thresholded.by_genes.txt") %>%
  column_to_rownames("Gene Symbol") %>%
  select_if(grepl("TARGET", names(.))) %>%
  distinct()#remove duplicate


CNV_lesions <- readr::read_tsv("../GenePattern GISTIC/GISTIC Purity Ploidy Q 0.25/OS_Purity_Ploidy_q_default.all_lesions.conf_90.txt") %>%
  rownames_to_column() %>%
  filter(!grepl("CN", `Unique Name`)) %>%
  column_to_rownames("Unique Name") %>%
  select_if(grepl("TARGET", names(.)))

CNV_lesions1 <- readr::read_tsv("../GenePattern GISTIC/GISTIC Purity Ploidy Q 0.25/OS_Purity_Ploidy_q_default.all_lesions.conf_90.txt") %>%
  rownames_to_column() %>%
  filter(grepl("CN", `Unique Name`)) %>%
  column_to_rownames("Unique Name") %>%
  select_if(grepl("TARGET", names(.)))

library(NMF)
estim.r <- nmf(d, 2:6,  nrun=500, seed=123, .opt='p6')
#save.image("Osteosarcoma_NMF.RData")
plot(estim.r)
###############
pdf("NMFclusters.pdf", width = 10, height = 12)
consensusmap(estim.r)
dev.off()

estim.r <- nmf(CNV_lesions, 2:6,  nrun=500, seed=123, .opt='p6')
pdf("NMFclusters Concensusmap.pdf", width = 10, height = 12)
#plot(estim.r)
consensusmap(estim.r)
dev.off()

pdf("NMF Estimate Plot.pdf", width = 10, height = 12)
plot(estim.r)
dev.off()

estim.r <- nmf(CNV_lesions, 2,  nrun=500, seed=123, .opt='p6')
predict(estim.r)
table(predict(estim.r))
pdf("NMFclusters Concensusmap Two Cluster.pdf", width = 10, height = 12)
#plot(estim.r)
consensusmap(estim.r)
dev.off()


### Concensus clusterPlus is not giving good results
CNV_lesions1 <- readr::read_tsv("../GenePattern GISTIC/GISTIC Purity Ploidy Q 0.25/OS_Purity_Ploidy_q_default.all_lesions.conf_90.txt") %>%
  rownames_to_column() %>%
  filter(!grepl("CN", `Unique Name`)) %>%
  column_to_rownames("Unique Name") %>%
  select_if(grepl("TARGET", names(.)))
dir.create("12348"); title="12348"
results.coding = ConsensusClusterPlus(as.matrix(CNV_lesions),maxK=10,reps=1000,pItem=1,pFeature=0.9, title=title,
                                      clusterAlg="hc",distance="euclidean", seed=12345,plot="png")
##############################################################

save.image("TARGET-OS-Cluatering.RData")
###################################################################
###################################################################