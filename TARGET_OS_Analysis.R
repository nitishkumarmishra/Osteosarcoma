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
library(matrixStats)
######################################################
setwd("D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/GDC Osteosarcoma")
#dir.create(file.path("TARGET-OS"), showWarnings = FALSE)

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
#matrixStats::rowVars(as.matrix(OS.HTSeq.Counts))

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
##############################################################################
## This part is to generate gencode.v22.gtf. I already save an loaded by using previous line
gtf <- rtracklayer::import('gencode.v22.annotation.gtf.gz')
gencode.v22.gtf <- as.data.frame(gtf)

gencode.v22.gtf.TCGA.selected <- gencode.v22.gtf %>%
  filter(type=="gene") %>%
  rename(Chr=seqnames, Type=type,ENSG=gene_id, GeneType=gene_type, Symbol=gene_name, Status=gene_status) %>%
  select(Chr, Type, GeneType, ENSG, Symbol)

saveRDS(gencode.v22.gtf.TCGA.selected, file = "gencode.v22.annotation.TCGA.rds")
#gencode.v34.annotation <- readRDS("D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/GDC Osteosarcoma/gencode.v34.annotation.rds")
count <- OS.HTSeq.FPKM
count <- count[rowSums(count==0)< ncol(count)*0.20,] ## Remove all gene which have 25% zero's
keep <- rowSums(edgeR::cpm(count)>1) >= ncol(count)*0.10 #### Select only genes which have have CPM > 1 for >=50% samples
count <- count[keep,]
count <- log2(OS.HTSeq.Counts[rownames(count),]+1)
############### gene median center normalization
############### https://www.sciencedirect.com/science/article/pii/S1535610817300016?via%3Dihub
#https://www.cell.com/cms/10.1016/j.ccell.2017.01.001/attachment/d7e976b2-aac3-4140-a993-5d212979a499/mmc1.pdf
# count <- log2(OS.HTSeq.FPKM+1)
#count <- sweep(count,1, apply(count,1,median,na.rm=T))
count <- count-apply(count, 1, median) # sweep and this give same value
#count <- t(scale(t(count), scale = FALSE))## both are same
# EnsembleID_Full <- rownames(count)
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
##### On Median Center Normalizes FPKM; MAD based #######
## https://bioconductor.org/packages/release/bioc/vignettes/ConsensusClusterPlus/inst/doc/ConsensusClusterPlus.pdf
library(ConsensusClusterPlus)

d.FPKM.coding <- as.matrix(dataset.protein.coding[,1:88])

dir.create("12345"); title="12345"
results.coding = ConsensusClusterPlus(d.FPKM.coding,maxK=10,reps=1000,pItem=1,pFeature=0.9, title=title,
                               clusterAlg="hc",distance="pearson", seed=12345,plot="png")

d.FPKM.lincRNA <- as.matrix(dataset.lincRNA[,1:88])
dir.create("12346"); title="12346"
results.lincRNA = ConsensusClusterPlus(d.FPKM.lincRNA,maxK=10,reps=1000,pItem=1,pFeature=0.9, title=title,
                               clusterAlg="hc",distance="pearson", seed=12345,plot="png")

################# NMF Clustering ########################
#########################################################
library(NMF)
d <- OS.HTSeq.FPKM
d <- d[rowSums(d==0)< ncol(d)*0.20,] ## Remove all gene which have 25% zero's
keep <- rowSums(edgeR::cpm(d)>1) >= ncol(d)*0.10 #### Select only genes which have have CPM > 1 for >=50% samples
d <- d[keep,]
d <- log2(d+1)
EnsembleID <- rownames(d)
d <- d %>%
  mutate(Mean= rowMeans(.), Median=rowMedians(as.matrix(.)) ,Var=rowVars(as.matrix(.)), 
       stdev=rowSds(as.matrix(.)), MAD=rowMads(as.matrix(.), method = "median")) %>%
  mutate(ENSG=EnsembleID)

Num_Gene=3000
d.protein.coding <- d %>%
  inner_join(gencode.v22.gtf.TCGA.selected, by = "ENSG") %>%
  filter(GeneType %in% c("protein_coding")) %>%
  #filter(GeneType %in% c("protein_coding"))%>% 
  slice_max(Var, n=Num_Gene)

Num_Gene=500
d.lincRNA <- d %>%
  inner_join(gencode.v22.gtf.TCGA.selected, by = "ENSG") %>%
  filter(GeneType %in% c("lincRNA")) %>%
  #filter(GeneType %in% c("protein_coding"))%>% 
  slice_max(MAD, n=Num_Gene)

d.protein.coding <- as.matrix(d.protein.coding[,1:88])
estim.r.protein <- nmf(d.protein.coding, 2:6,  nrun=500, seed=123, .opt='p6')
pdf("NMFclusters HTSeq FPKM Rank Survey.pdf", width = 10, height = 12)
plot(estim.r.protein)
dev.off()

pdf("NMFclusters HTSeq FPKM.pdf", width = 10, height = 12)
consensusmap(estim.r.protein)
dev.off()

si.1 <- silhouette(nmf(d.protein.coding, 4, nrun = 500))
summary(si.1) ## To check silhouette width
plot(si.1)
plot(si.1, col = c("red", "blue", "green", "yellow"))

d.lincRNA <- as.matrix(d.lincRNA[,1:88])
estim.r.lincRNA <- nmf(d.lincRNA, 2:6,  nrun=500, seed=123, .opt='p6')
pdf("NMFclusters HTSeq lincRNA FPKM Rank Survey.pdf", width = 10, height = 12)
plot(estim.r.lincRNA)
dev.off()

pdf("NMFclusters HTSeq lincRNA FPKM.pdf", width = 10, height = 12)
consensusmap(estim.r.lincRNA)
dev.off()

si.2 <- silhouette(nmf(d.lincRNA, 4, nrun = 1000))
summary(si.2) ## To check silhouette width
plot(si.2)
plot(si.2, col = c("red", "blue", "green", "yellow"))

#########################################################
############# Somatic mutation analysis #################
############# For all samples ###########################
library(maftools)
Osteo <- read.maf(maf = "Combined_All_MuTect2.aliquot.maf", clinicalData = NULL)
Osteo.SomaticSniper<- read.maf(maf = "Combined_All_SomaticSniper.aliquot.maf", clinicalData = NULL)
Osteo.MuSE <- read.maf(maf = "Combined_All_MuSE.aliquot.maf", clinicalData = NULL)
Osteo.VarScan2 <- read.maf(maf = "Combined_All_VarScan2.aliquot.maf", clinicalData = NULL)
Osteo.Pindel <- read.maf(maf = "Combined_All_Pindel.aliquot.maf", clinicalData = NULL)
Osteo.ensemble_raw <- read.maf(maf = "Combined_All_aliquot_ensemble_raw.maf", clinicalData = NULL)
Osteo.ensemble_masked <- read.maf(maf = "Combined_All_aliquot_ensemble_masked.maf", clinicalData = NULL)
getSampleSummary(Osteo)
getGeneSummary(Osteo)
plotmafSummary(maf = Osteo, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
oncoplot(Osteo, top = 10)
oncoplot(maf = Osteo, draw_titv = TRUE)
lollipopPlot(maf = Osteo, gene = 'PABPC3', AACol = 'HGVSp_Short', showMutationRate = TRUE)
Osteo.titv = titv(maf = Osteo, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = Osteo.titv)
#OncogenicPathways function checks for enrichment of known Oncogenic Signaling Pathways in TCGA cohorts
#Sanchez-Vega F,  et al. Oncogenic Signaling Pathways in The Cancer Genome Atlas. Cell (2018) 173: 321-337 e310
OncogenicPathways(maf = Osteo)
# drugInteractions function checks for drug-gene interactions and gene druggability information compiled from Drug Gene Interaction database.
dgi = drugInteractions(maf = Osteo, fontSize = .5)

#########################################################
############# Somatic mutation analysis #################
######### Common in Mutation (MAF) and RNAseq ###########
## Make clinical and MAF file for samples which are common in MAF and RNAseq
Trial <- Osteo
Trial@clinical.data$Tumor_Sample_Barcode <- substr(Trial@clinical.data$Tumor_Sample_Barcode, 1, 16)
Trial@data$Tumor_Sample_Barcode <- substr(Trial@data$Tumor_Sample_Barcode, 1, 16)
Trial@data$Matched_Norm_Sample_Barcode <- substr(Trial@data$Matched_Norm_Sample_Barcode, 1, 16)
Trial@variants.per.sample$Tumor_Sample_Barcode <- substr(Trial@variants.per.sample$Tumor_Sample_Barcode, 1, 16)
Trial@variant.type.summary$Tumor_Sample_Barcode <- substr(Trial@variant.type.summary$Tumor_Sample_Barcode, 1, 16)
Trial@maf.silent$Tumor_Sample_Barcode <- substr(Trial@maf.silent$Tumor_Sample_Barcode, 1, 16)
# Osteo.RNAseq <- subsetMaf(maf = Trial, tsb = substr(colnames(OS.HTSeq.Counts), 1, 16), isTCGA = FALSE, clinQuery=TRUE)
# ttt <- getClinicalData(x = Trial)
# Osteo.RNAseq@clinical.data <- ttt[ttt$Tumor_Sample_Barcode %in% substr(colnames(OS.HTSeq.Counts), 1, 16),]
#Osteo.RNAseq <- subsetMaf(maf = Osteo, tsb = substr(colnames(OS.HTSeq.Counts), 1, 16), isTCGA = FALSE, clinQuery=TRUE)
ttt <- Osteo@clinical.data
ttt$TCGA_short <- substr(ttt$Tumor_Sample_Barcode, 1, 16)
TARGET_ID <- intersect(ttt$TCGA_short, substr(colnames(OS.HTSeq.Counts), 1, 16))
OS.Clinical.RNAseq <- OS.Clinical[OS.Clinical$Tumor_Sample_Barcode %in%TARGET_ID,]

#ttt$Extra <- substr(ttt$Tumor_Sample_Barcode, 1, 16)
colnames(ttt) <- c("Extra","Tumor_Sample_Barcode")
tttt <- merge(OS.Clinical.RNAseq, ttt, by="Tumor_Sample_Barcode")
tttt <- tttt[!duplicated(tttt$Tumor_Sample_Barcode),] ## Few patients have two file.
tttt$Tumor_Sample_Barcode <- NULL
tttt$Tumor_Sample_Barcode <- tttt$Extra; tttt$Extra <- NULL
write.csv(tttt$Tumor_Sample_Barcode, file = "RNAseqID.txt", row.names = FALSE)

# Remove header and " from the file, copy on rhino 
# dos2unix RNAseqID.txt
# head -n 8 Combined_All_MuTect2.aliquot.maf >RNAseq_Header
# grep -Fwf RNAseqID.txt Combined_All_MuTect2.aliquot.maf >tmp
# cat RNAseq_Header tmp >Combined_RNAseq_MuTect2.aliquot.maf
# grep -Fwf RNAseqID.txt Combined_All_SomaticSniper.aliquot.maf >tmp
# cat RNAseq_Header tmp >Combined_RNAseq_SomaticSniper.aliquot.maf
# grep -Fwf RNAseqID.txt Combined_All_VarScan2.aliquot.maf >tmp
# cat RNAseq_Header tmp >Combined_RNAseq_VarScan2.aliquot.maf

#########################################################
Osteo.RNAseq <- read.maf(maf = "Combined_RNAseq_MuTect2.aliquot.maf", clinicalData = tttt)
getSampleSummary(Osteo.RNAseq)
getGeneSummary(Osteo.RNAseq)
plotmafSummary(maf = Osteo.RNAseq, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
oncoplot(Osteo.RNAseq, top = 10)
oncoplot(maf = Osteo.RNAseq, draw_titv = TRUE)
lollipopPlot(maf = Osteo.RNAseq, gene = 'PABPC3', AACol = 'HGVSp_Short', showMutationRate = TRUE)
Osteo.titv = titv(maf = Osteo.RNAseq, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = Osteo.titv)
#OncogenicPathways function checks for enrichment of known Oncogenic Signaling Pathways in TCGA cohorts
#Sanchez-Vega F,  et al. Oncogenic Signaling Pathways in The Cancer Genome Atlas. Cell (2018) 173: 321-337 e310
OncogenicPathways(maf = Osteo)
# drugInteractions function checks for drug-gene interactions and gene druggability information compiled from Drug Gene Interaction database.
dgi = drugInteractions(maf = Osteo, fontSize = .5)

#########################################################
############# Clustering mutation analysis###############
concensusName <- names(results.coding[[2]]$consensusClass)
concensusClusterPlus.Cluster1 <- concensusName[results.coding[[2]]$consensusClass==1]
concensusClusterPlus.Cluster2 <- concensusName[results.coding[[2]]$consensusClass==2]
concensusClusterPlus.Cluster1 <- substr(concensusClusterPlus.Cluster1, 1, 16)
concensusClusterPlus.Cluster2 <- substr(concensusClusterPlus.Cluster2, 1, 16)
write.table(concensusClusterPlus.Cluster1, file = "concensusClusterPlus.Cluster1.txt", row.names = FALSE, col.names = FALSE)
write.table(concensusClusterPlus.Cluster2, file = "concensusClusterPlus.Cluster2.txt", row.names = FALSE, col.names = FALSE)
#dos2unix concensusClusterPlus.Cluster1.txt concensusClusterPlus.Cluster2.txt 
#grep -Fwf concensusClusterPlus.Cluster1.txt Combined_RNAseq_MuTect2.aliquot.maf >tmp
#cat RNAseq_Header tmp >Combined_RNAseq_MuTect2.concensusClusterPlus.Cluster1.maf
#grep -Fwf concensusClusterPlus.Cluster2.txt Combined_RNAseq_MuTect2.aliquot.maf >tmp
#cat RNAseq_Header tmp >Combined_RNAseq_MuTect2.concensusClusterPlus.Cluster2.maf

concensusClusterPlus.Cluster1.maf = read.maf(maf = "Combined_RNAseq_MuTect2.concensusClusterPlus.Cluster1.maf", clinicalData = NULL)
concensusClusterPlus.Cluster2.maf = read.maf(maf = "Combined_RNAseq_MuTect2.concensusClusterPlus.Cluster2.maf", clinicalData = NULL)

#Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
Class1.vs.Class2 <- mafCompare(m1 = concensusClusterPlus.Cluster1.maf, m2 = concensusClusterPlus.Cluster2.maf, m1Name = 'Class1', m2Name = 'Class2', minMut = 5)
print(Class1.vs.Class2)
forestPlot(mafCompareRes = Class1.vs.Class2, pVal = 0.025, color = c('royalblue', 'maroon'), geneFontSize = 0.8)

genes = c("MUC5AC", "INO80D","ZNF880", "CELF3", "SACS","SIPA1L1",  "MACC1", "PRG2")
coOncoplot(m1 = concensusClusterPlus.Cluster1.maf , m2 = concensusClusterPlus.Cluster2.maf , m1Name = 'Class1-OS', m2Name = 'Class2-OS', genes = genes, removeNonMutated = TRUE)

## Mutually exclusive or co-occurring set of genes can be detected using somaticInteractions function, 
## Which performs pair-wise Fisher's Exact test to detect such significant pair of genes.
## https://github.com/jhrcook/wext
## http://compbio.cs.brown.edu/projects/wext/
somaticInteractions(maf = Osteo.RNAseq, top = 25, pvalue = c(0.05, 0.1))

# Detecting cancer driver genes based on positional clustering
#maftools has a function oncodrive which identifies cancer genes (driver) from a given MAF. oncodrive is a based on algorithm oncodriveCLUST
Osteo.RNAseq.sig = oncodrive(maf = Osteo.RNAseq, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
head(Osteo.RNAseq.sig)
plotOncodrive(res = Osteo.RNAseq.sig, fdrCutOff = 0.5, useFraction = FALSE)

#Adding and summarizing pfam domains
#maftools comes with the function pfamDomains, which adds pfam domain information to the amino acid changes
Osteo.RNAseq.pfam = pfamDomains(maf = Osteo.RNAseq, AACol = 'HGVSp_Short', top = 10)
#Protein summary (Printing first 7 columns for display convenience)
Osteo.RNAseq.pfam$proteinSummary[,1:7, with = FALSE]
#Domain summary (Printing first 3 columns for display convenience)

Osteo.RNAseq.pfam$domainSummary[,1:3, with = FALSE]
########################################################
############# GenVisR waierfall plot ###################
library(GenVisR)
set.seed(426)
mafObject <- MutationAnnotationFormat("Combined_RNAseq_MuTect2.aliquot.maf", version = 2.4)
drawPlot(Waterfall(mafObject, recurrence = 0.8))
drawPlot(Waterfall(mafObject, recurrence = 0.8, geneMax = 20))

############# DEG analysis and GSEA #####################
#########################################################
Expression <- OS.HTSeq.Counts
Expression <- Expression[rowSums(Expression==0)< ncol(Expression)*0.20,] ## Remove all gene which have 25% zero's
keep <- rowSums(edgeR::cpm(Expression)>1) >= ncol(Expression)*0.10 #### Select only genes which have have CPM > 1 for >=50% samples
Expression <- Expression[keep,]

EnsembleID <- rownames(Expression)
Expression$ENSG <- EnsembleID

Expression.protein.coding <- Expression %>%
  inner_join(gencode.v22.gtf.TCGA.selected, by = "ENSG") %>%
  filter(GeneType %in% c("protein_coding")) 
Expression.protein.coding <- Expression.protein.coding[!duplicated(Expression.protein.coding$Symbol),]
rownames(Expression.protein.coding) <- Expression.protein.coding$Symbol

Expression.protein.coding <- Expression.protein.coding[,grepl("TARGET", colnames(Expression.protein.coding))]
colnames(Expression.protein.coding) <- substr(colnames(Expression.protein.coding), 1, 16)
Expression.protein.coding.Group1 <- Expression.protein.coding[,concensusClusterPlus.Cluster1]
Expression.protein.coding.Group2 <- Expression.protein.coding[,concensusClusterPlus.Cluster2]

Expression.Group <- cbind(Expression.protein.coding.Group1, Expression.protein.coding.Group2)


library(DESeq2)
condition <- factor(c(rep("Group1", 45),rep("Group2", 43)), levels = c("Group1", "Group2"))
countdata <- Expression.Group
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
#dds
dds <- DESeq(dds)
res <- results(dds)
table(res$padj<0.05)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=FALSE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
resdata.p.0.05 <- resdata[resdata$padj <= 0.05,]

library(clusterProfiler)
library(org.Hs.eg.db)
detach("package:dplyr", unload=TRUE)

gene <- resdata.p.0.05$Gene
gene <- sort(gene, decreasing = TRUE)
gmtfile <- system.file("extdata", "c2.cp.kegg.v7.1.symbols.gmt", package="clusterProfiler")
c2 <- read.gmt(gmtfile)
egmt <- enricher(gene, TERM2GENE=c2, pvalueCutoff = 0.01, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2)
summary(egmt)

dotplot(egmt, showCategory=10)
boxplot(egmt, showCategory=10)


ranks <- resdata.p.0.05$log2FoldChange
names(ranks) <- resdata.p.0.05$Gene
head(ranks)
ranks <- sort(ranks, decreasing = T)
barplot(sort(ranks, decreasing = T))


library(fgsea)
library(data.table)
library(ggplot2)
C2 <- qusage::read.gmt(gmtfile)
fgseaRes <- fgsea(pathways = C2, 
                  stats    = ranks,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)
head(fgseaRes[order(pval), ])

plotEnrichment(C2[["KEGG_OXIDATIVE_PHOSPHORYLATION"]],ranks, ticksSize = 1) +
  labs(title="Oxydative Phosphorylation") +
  geom_point(color = "brown", size = 1.5) +
  geom_line(size = 1, color = "green4") 
   
plotEnrichment(C2[["KEGG_CALCIUM_SIGNALING_PATHWAY"]],ranks, ticksSize = 1) +
    labs(title="Calcium Signaling Pathway") +
    geom_point(color = "brown", size = 1.5) +
    geom_line(size = 1, color = "green4")

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(C2[topPathways], ranks, fgseaRes, gseaParam=0.5)

collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], C2, ranks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
plotGseaTable(C2[mainPathways], ranks, fgseaRes, gseaParam = 0.5)

#########################################################
############### Prepare for MutSigCV ####################
OS.mutsig.corrected = prepareMutSig(maf = Osteo.RNAseq)
Cluster1.mutsig.corrected <- prepareMutSig(maf=concensusClusterPlus.Cluster1.maf)
Cluster2.mutsig.corrected <- prepareMutSig(maf=concensusClusterPlus.Cluster2.maf)

write.table(OS.mutsig.corrected, file = "OS.RNAseq.mutsig.maf", quote=FALSE, sep='\t', row.names = FALSE)
write.table(Cluster1.mutsig.corrected, file = "Cluster1.mutsig.corrected.maf", quote=FALSE, sep='\t', col.names = NA)
write.table(Cluster2.mutsig.corrected, file = "Cluster2.mutsig.corrected.maf", quote=FALSE, sep='\t', col.names = NA)

#########################################################
################### TARGET-OS CNV #######################
library(TargetOsteoAnalysis)

variants = variant_calls()
#OS.GISTIC <- target_gistic_se()
gistic_file = system.file(package='TargetOsteoAnalysis','extdata/all_data_by_genes.txt.bz2')
OS.GISTIC = read_tsv(gistic_file)
OS.Purity = target_purity()
#x = target_cnv_calls()
# pureCN of TARGET OS
fname = system.file("extdata/summary.pureCN.cns.calls.tsv.gz", package = "TargetOsteoAnalysis")
OS.CNVkit.Purity = read_tsv(fname, col_names = TRUE)

# data(variation.hg18.v10.nov.2010)
# gbm.cn=CNregions(seg=gbm.seg,epsilon=0,adaptive=FALSE,rmCNV=TRUE,
#                    cnv=variation.hg18.v10.nov.2010[,3:5],
#                    frac.overlap=0.5, rmSmallseg=TRUE,nProbes=5)
## Chromosome arm level length
library(data.table)
Cytoband <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz",
col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
Cytoband <- Cytoband[ , .(length = sum(chromEnd - chromStart)), by = .(chrom, arm = substring(name, 1, 1)) ]
Cytoband <- Cytoband[grep("p|q", Cytoband$arm),]
Cytoband <- Cytoband[gtools::mixedorder(Cytoband$chrom),]

Gene.Cytoband <- OS.GISTIC %>% select(c("Gene Symbol", "Cytoband"))
Gene.Cytoband$GeneSymbol <- sapply(strsplit(Gene.Cytoband$`Gene Symbol`, "\\|"),'[[',1)
Gene.Cytoband$`Gene Symbol` <- NULL

OS.CNVkit.Purity$GeneSymbol <- sapply(strsplit(OS.CNVkit.Purity$gene, "\\,"),'[[',1)
CNV.Position <- OS.CNVkit.Purity %>%
  select("chromosome", "start", "end")

genes.cytoband <- fread("gene_bandall.txt")
genes.cytoband <- genes.cytoband %>% 
  select(c("V2", "V6", "V5"))
genes.cytoband$V5 <- gsub("\\-|\\+| ", "", genes.cytoband$V5)
colnames(genes.cytoband) <- c("GeneSymbol", "Arm", "Coordinate")
genes.cytoband <- unique(genes.cytoband)

Cytoband <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz",
col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
Cytoband <- Cytoband[grep("p|q", Cytoband$name),]
Cytoband <- Cytoband[gtools::mixedorder(Cytoband$chrom),]
write.csv(Cytoband, file = "Cytoband.bed")

## File generated in CBSB by using bedtools (details in Evernote)
CNV.Purity <- read_tsv("CNV_Cytoband.txt", col_names = c("Pchr", "Pstart", "Pend", "strand", 
                                                         "arm",  "chromosome", "start", "end")) %>%
  select(chromosome, start, end, arm, Pstart, Pend) %>%
  mutate(ArmLength = Pend-Pstart)

CNV.Purity <- unique(CNV.Purity)
CNV.Purity.1 <- inner_join(OS.CNVkit.Purity,CNV.Purity, by=c("chromosome", "start", "end")) %>%
  mutate(ArmRatio=CNregionLength/ArmLength)

CNV.Purity.2 <- CNV.Purity.1[CNV.Purity.1$ArmRatio >= 0.5,]
CNV.Purity.2 <- CNV.Purity.2 %>% 
  unite("CHR", c(chromosome,arm), sep = ".", remove = FALSE)

CNV.Purity.3 <- CNV.Purity.2 %>%
  select(CHR, dataName, log2) %>%
  group_by(CHR, dataName) %>%
  summarise(mean=mean(log2)) %>%
  #mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = dataName, values_from = mean, values_fill = 0) 
  #select(-row)


Segment <- readr::read_tsv("../GenePattern GISTIC/TARGET_OS_L3_Segmentation.txt", col_names = TRUE) %>%
  mutate(Chromosome = gsub("chrX", "chr23", Chromosome)) %>%
  mutate(Chromosome = gsub("chrY", "chr24", Chromosome)) %>%
  filter(Chromosome != "chrMT") %>% 
  mutate(Chromosome = gsub("chr", "", Chromosome)) %>%
  mutate(Chromosome = as.numeric(Chromosome))

#dir = "D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/GenePattern GISTIC/"
#readr::write_tsv(Segment, file.path(dir, "Segment.tsv")) ## Use it for GISTIC (GenePattern)
#TCGA.CNV.bypass <- readr::read_tsv("../GenePattern GISTIC/CNV.hg19.bypos.111213.txt", col_names = TRUE) %>%
 # rename(chr=Chromosome, start=Start, end=End) %>%
 # mutate(Chromosome=(paste0("chr", as.character(Chromosome))))
#TCGA.CNV.bypass <- as.matrix(TCGA.CNV.bypass)
TCGA.CNV.bypass <- read.csv("../GenePattern GISTIC/CNV.hg19.bypos.111213.txt", header = TRUE, sep = "\t")
TCGA.CNV.bypass$Chromosome <- paste0("chr", TCGA.CNV.bypass$Chromosome)

library(iClusterPlus)
library(GenomicRanges)
library(gplots)
library(lattice)
library(cluster)
seg <- as.data.frame(Segment)
OS.cn=CNregions(seg=seg,epsilon=0.02,adaptive=FALSE,rmCNV=TRUE,
                 cnv=TCGA.CNV.bypass[,2:4],
                 frac.overlap=0.5, rmSmallseg=TRUE,nProbes=5) # Total 6019 regions

#OS.cn=CNregions(seg=seg,adaptive=TRUE,rmCNV=TRUE, cnv=TCGA.CNV.bypass[,2:4],
#                frac.overlap=0.5, rmSmallseg=TRUE,nProbes=10) ## Decreased number of regions 2122

Osteo.MuTect.GISTIC <- read.maf(maf = "Combined_All_MuTect2.aliquot.maf", clinicalData = NULL,
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

######### https://www.nature.com/articles/ng.2760.pdf
##alpha=purity; tau=ploidy
# Rx = (alpha*Qx+2(1-alpha))/D
# D = (alpha*tau) + 2(1-alpha)
# Qx = D-Rx/alpha -2(1-alpha)/alpha
# Rx.1 = qx/tau
# Rx.1 = (Rx/alpha) - 2*(1-alpha)/(alpha*tau)

Segment <- Segment %>%
  mutate(Source = substr(Segment$SampleName, 11, 16))

Segment.Purity <- left_join(Segment, OS.Purity, by=c("Source")) %>%
  select(!c(Source, dataset)) %>%
  mutate("Segment Purity" = if_else(!is.na(Purity),((`Segment Mean`/Purity) - 2*(1-Purity)/(Purity*Ploidy)),`Segment Mean`))

Segment.Purity.selected <- Segment.Purity %>%
  select(c(SampleName, Chromosome, Start, End, `Num Probes`, `Segment Purity`)) %>%
  rename("Segment Mean"="Segment Purity")

write.table(Segment.Purity.selected, file = "../GenePattern GISTIC/TARGET_OS_L3_Segmentation.Purity.txt", col.names = TRUE, row.names = FALSE, sep = "\t")
## Replace " with NA in notepad++

Segment.Purity.0.5 <- Segment.Purity %>% group_by(SampleName) %>% filter(Purity >= 0.5) %>%
  select(c(SampleName, Chromosome, Start, End, `Num Probes`, `Segment Purity`)) %>%
  dplyr::rename("Segment Mean"="Segment Purity")

  
write.table(Segment.Purity.0.5, file = "../GenePattern GISTIC/TARGET_OS_L3_Segmentation.Purity.0.5.txt", col.names = TRUE, row.names = FALSE, sep = "\t")
## Replace " with NA in notepad++

seg.purity <- as.data.frame(Segment.Purity.selected)
OS.cn.purity <- CNregions(seg=seg.purity,epsilon=0.05,adaptive=TRUE,rmCNV=TRUE,
                cnv=TCGA.CNV.bypass[,2:4],
                frac.overlap=0.5, rmSmallseg=TRUE,nProbes=5) # Total 2022 regions

OS.cn.purity.1 <- CNregions(seg=seg.purity,epsilon=0.04,adaptive=FALSE,rmCNV=TRUE,
                          cnv=TCGA.CNV.bypass[,2:4],
                          frac.overlap=0.6, rmSmallseg=TRUE,nProbes=5) # Total 5996 regions
#########################################
df <- CNV.Purity.3 %>%
  column_to_rownames("CHR")
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms

hc <- agnes(t(df), method = "ward")
#pltree(hc, cex = 0.6, hang = -1, main = "Dendrogram of Chromosome Arms") 
d <- dist(t(df), method = "euclidean")
hc1 <- hclust(d, method = "ward.D" )
sub_grp <- cutree(hc1, k = 3)

plot(hc1, cex = 0.6)
rect.hclust(hc1, k = 3, border = 2:5)

## change k in cutree to check best shilhoutte plot 
sil_cl <- silhouette(cutree(hc1, k=2, h=25) ,d, title=title(main = 'Good'))
rownames(sil_cl) <- rownames(d)
plot(sil_cl, col = c("red", "blue"))
d1 <- cophenetic(hc1)
cor(d, d1) ## This is cophenetic correlation coefficients

##miRNA expression of TARGET-OS
## With normalization
genedat = "ftp://caftpd.nci.nih.gov/pub/OCG-DCC/TARGET/OS/miRNA_pcr/L2/TARGET_OS_miRNA_level2_dCt.txt"
OS.miRNA = readr::read_tsv(genedat)
## without normalization miRNA 
genedat = "https://target-data.nci.nih.gov/Public/OS/miRNA_pcr/L2/TARGET_OS_miRNA_level2_Ct.txt"
OS.miRNA1 = readr::read_tsv(genedat)


d1 <- OS.miRNA1
d1 <- d1[!duplicated(d1$miRNA),]
miRNA <- d1$miRNA
d1 <- d1 %>%
  tibble::column_to_rownames("miRNA")

d1 <- d1 %>%
  mutate(Mean= rowMeans(.), Median=rowMedians(as.matrix(.)) ,Var=rowVars(as.matrix(.)), 
         SD=rowSds(as.matrix(.)), MAD=rowMads(as.matrix(.), method = "median")) %>%
  mutate(miRNA=miRNA)

Num_Gene=nrow(d1)*0.25
d1.miRNA <- d1 %>%
    slice_max(Var, n=Num_Gene)

d1.NMF <- as.matrix(d1.miRNA[,1:89])
estim.r.miRNA <- NMF::nmf(d1.NMF, 2:10,  nrun=500, seed=123, .opt='p6')
pdf("NMFclusters miRNA Rank Survey.pdf", width = 10, height = 12)
plot(estim.r.miRNA)
dev.off()

pdf("NMFclusters miRNA.pdf", width = 10, height = 12)
consensusmap(estim.r.miRNA)
dev.off()


si.1 <- silhouette(nmf(d1.NMF, 3, nrun = 500))
summary(si.1) ## To check silhouette width
plot(si.1)
plot(si.1, col = c("red", "blue", "green"))
#plot(si.1, col = c("red", "blue", "green", "yellow"))
#########################################################
#################### CPDB GSEA data #####################
# The gene set collection from CPDB can be formatted by: 
# "CPDB_pathways_genes.tab" is available in the download section in http://consensuspathdb.org/
CPDB = readLines("D:/OneDrive - University of Nebraska Medical Center/MSigDB/CPDB_pathways_genes.tab", warn = FALSE)
CPDB = lapply(CPDB, function(x) strsplit(x, "\t")[[1]])
names(CPDB) = sapply(CPDB, function(x) paste(x[3], x[2], x[1], sep=":") )
CPDB = lapply(CPDB, function(x) x[-c(1:3)])
CPDB = CPDB[-length(CPDB)]
CPDB = lapply(CPDB, function(x) strsplit(x, ",")[[1]])
CPDB$`source:external_id:pathway` <- NULL

#########################################################
################### DNA Methylation #####################
OS.Meth <- readRDS("../OsteoMeth.rds")
hm450.tss1500.Masked.Uniq <- readRDS("../hm450.tss1500.Masked.Uniq.rda")

#########################################################
################### Save workspace ######################
rm(query, gtf, EnsembleID, keep, EnsembleID_Full, gistic_file, fname, genedat)
## First time I save rworkspace of read counts and clinical data in GDC-Harmonized folder in eMAP
#save.image("D:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/GDC-Harmonized/TARGET-OS/OS.RData")
#save.image("Concensus_Clustering.RData")
#save.image("TARGET-OS.RData")
rm(list = ls())
##########################################################
save(OS.Meth, hm450.tss1500.Masked.Uniq, OS.FPKM.UQ, OS.HTSeq.FPKM, OS.Clinical.TARGET, 
     gencode.v22.gtf.TCGA.selected, seg.purity, OS.miRNA, OS.miRNA1, d1.miRNA, CNV.Purity.3, MuTect.matrix, 
     MuTect.matrix.Binary2, OS.cn.purity, OS.cn.purity.1, file = "Osteo_iCluster.RData")
## Find common samples across different types of the data

Common.sample <- Reduce(intersect, list(substr(OS.Clinical.TARGET$`TARGET USI`, 11, 16), #Clinical
                       substr(colnames(d1.miRNA), 11, 16), # miRNA
                       substr(colnames(CNV.Purity.3), 1, 6), #CNV
                       substr(rownames(OS.cn.purity.1), 11,16), #CNV
                       substr(rownames(MuTect.matrix.Binary2), 11,16),#SNV
                       substr(colnames(OS.FPKM.UQ), 11,16),#RNAseq
                       substr(colnames(OS.Meth), 11,16))) #Methylation

Common.sample.purity.0.5 <- Reduce(intersect, list(substr(OS.Clinical.TARGET$`TARGET USI`, 11, 16), #Clinical
                       substr(colnames(d1.miRNA), 11, 16), # miRNA
                       substr(colnames(CNV.Purity.3), 1, 6), #CNV
                       substr(rownames(OS.cn.purity.1), 11,16), #CNV
                       substr(rownames(MuTect.matrix.Binary2), 11,16),#SNV
                       substr(colnames(OS.FPKM.UQ), 11,16),#RNAseq
                       substr(colnames(OS.Meth), 11,16)),#Methylation
                       substr(Segment.Purity.0.5$SampleName, 11, 16)) # Purity >0.5

## CNV based clustering https://www.sciencedirect.com/science/article/pii/S1535610817302957
## gene expression normalization based on purity
## https://bmcmedgenomics.biomedcentral.com/articles/10.1186/1755-8794-4-54

#ifelse(grepl("TARGET", colnames(d)), substr(colnames(d), 11,16), colnames(d))
OS.Meth <- OS.Meth[,grepl("TARGET", colnames(OS.Meth))]
colnames(OS.Meth) <- substr(colnames(OS.Meth), 11,16)


OS.Meth.Common <- OS.Meth[,Common.sample]
OS.Meth.Common <- t(OS.Meth.Common)

OS.Purity.1 <- OS.Purity %>% 
  arrange(desc(Purity)) %>% 
  distinct(Source, .keep_all = TRUE) %>%
  column_to_rownames(var = "Source")

#OS.Meth.Common %>% inner_join(OS.cn.purity.1, by="Source")
OS.Meth.Common <- merge(OS.Meth.Common, OS.Purity.1, by="row.names")
OS.Meth.Common <- OS.Meth.Common[, -grep("rs", colnames(OS.Meth.Common))]
OS.Meth.Common <- as.data.frame(OS.Meth.Common)
rownames(OS.Meth.Common) <- OS.Meth.Common$Row.names; OS.Meth.Common$Row.names <- NULL

meth <- OS.Meth.Common[, !colnames(OS.Meth.Common) %in% c("Purity", "Ploidy", "dataset")]
#meth <- as.matrix(meth)
purity <- OS.Meth.Common$Purity

set.seed(123)
ll <- mapply(function(x,y)cor.test(as.numeric(meth[,x]),purity, method = "spearman", alternative = "two.sided"),
             1:ncol(meth),
             SIMPLIFY=FALSE)

cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
spearman.p <- cbind(cor.value, p.value)
rownames(spearman.p) <- colnames(meth)
spearman.p <- data.frame(spearman.p)
spearman.p$BH <- p.adjust(spearman.p$p.value, method = "BH")
set.seed(123)
ll.perm <- mapply(function(x,y)jmuOutlier::perm.cor.test(as.numeric(meth[,x]),purity, method = "spearman", alternative = "two.sided", num.sim = 1000),
              1:ncol(meth),
              SIMPLIFY=FALSE)
 
perm.p.value <- sapply(ll.perm,'[[','p.value')
spearman.p <- cbind(spearman.p, perm.p.value)

purity.CpGs <- spearman.p[spearman.p$BH <= 0.01 & spearman.p$perm.p.value < 0.01 & abs(spearman.p$cor.value)>=0.25,]

meth.selected <- meth[,!colnames(meth) %in% rownames(purity.CpGs)]

keep <- colSums(meth.selected >0.2) >= nrow(meth.selected)*0.05 ## Select CpGs beta >0.2 for >10% samples
meth.selected <- meth.selected[,keep]


library(sesame)
sesame <- sesameDataGet('HM450.probeInfo')
SexChr <- sesame$mapped.probes.hg19[seqnames(sesame$mapped.probes.hg19)=="chrX" | seqnames(sesame$mapped.probes.hg19)=="chrY"]

library(magrittr)
meth.selected <- meth.selected %>%
  t %>%
  as.data.frame

ProbeID <- rownames(meth.selected)
meth.selected <- meth.selected %>%
  mutate(Mean= rowMeans(.), Median=rowMedians(as.matrix(.)) ,Var=rowVars(as.matrix(.)), 
         stdev=rowSds(as.matrix(.)), MAD=rowMads(as.matrix(.), method = "median")) %>%
  mutate(IlluminaID=ProbeID)

Num_Gene=3000
meth.cluster <- meth.selected %>%
  slice_max(Var, n=Num_Gene)


d.meth.cluster <- as.matrix(meth.cluster[,1:77])

dir.create("ConcensusMeth"); title="ConcensusMeth"
results.meth = ConsensusClusterPlus(d.meth.cluster,maxK=10,reps=1000,pItem=1,pFeature=0.9, title=title,
                                      clusterAlg="hc",distance="pearson", seed=12345,plot="png")


estim.r.meth <- NMF::nmf(d.meth.cluster, 2:10,  nrun=500, seed=123, .opt='p6')
pdf("NMFclusters Methylation Rank Survey.pdf", width = 10, height = 12)
plot(estim.r.meth)
dev.off()


si <- silhouette(nmf(d.meth.cluster, 4, nrun = 500))
summary(si) ## To check silhouette width
plot(si)
plot(si, col = c("red", "blue", "green", "purple"))
#############################################################
save.image("Concensus_Clustering.RData")