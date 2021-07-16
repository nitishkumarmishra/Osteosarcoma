setwd("D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma")
load("D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/OsteosarcomaData.RData")
load("D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/TARGET-GDC-OS-RNAseq.RData")
################################################################
################################################################
#common_sample <- intersect(substr(colnames(myNorm), 1, 20),substr(colnames(OS_HTSeq_FPKM_UQ), 1, 20))
Beta_Value <- myNorm
colnames(Beta_Value) <- substr(colnames(Beta_Value), 1, 20)

FPKM_UQ_Value <- OS_HTSeq_FPKM_UQ
colnames(FPKM_UQ_Value) <- substr(colnames(FPKM_UQ_Value), 1, 20)
rownames(FPKM_UQ_Value) <- gsub("\\.[1-9]", "", rownames(FPKM_UQ_Value))

common_sample <- intersect(colnames(Beta_Value), colnames(FPKM_UQ_Value))
Beta_Value.impute <- Beta_Value

Beta_Value.impute <- Beta_Value.impute[, common_sample]
Beta_Value.impute <- impute::impute.knn(Beta_Value.impute, k = 15, rowmax = 0.75, colmax = 0.75, rng.seed = 123)
Beta_Value.impute <- Beta_Value.impute$data

variance <- matrixStats::rowVars(Beta_Value.impute)
SD <- matrixStats::rowSds(Beta_Value.impute)
Beta_Value.impute <- cbind(Beta_Value.impute, variance, SD)
Beta_Value.impute <- as.data.frame(Beta_Value.impute)

Beta_Value.impute <- Beta_Value.impute[order(Beta_Value.impute$variance, decreasing = TRUE), ]

Beta_Value.impute.2000 <- Beta_Value.impute[Beta_Value.impute$SD >= 0.3,]

Beta_Value.impute.2000 <- Beta_Value.impute.2000[-grep("rs", rownames(Beta_Value.impute.2000)),]
Beta_Value.impute.2000$variance <- NULL; Beta_Value.impute.2000$SD <- NULL

d <- as.matrix(Beta_Value.impute.2000)
#########################################################
#########################################################
library(NMF)
estim.r <- nmf(d, 2:6,  nrun=500, seed=123, .opt='p6')
#save.image("Osteosarcoma_NMF.RData")
plot(estim.r)
###############
pdf("NMFclusters.pdf", width = 10, height = 12)
consensusmap(estim.r)
dev.off()
#####Select best K number, then calculate silhouette width ####
## We observed that consensusmap suggests two cluster have best result.
## Check the cophenetic coefficients at two cluser and three. We observed that coefficient value goes down at three cluster solution.
## So there is two possible subtype of the Osteosarcoma.
si <- silhouette(nmf(d, 3, nrun = 500))
summary(si) ## To check silhouette width
plot(si)
plot(si, col = c("red", "blue", "green"))

# si <- silhouette(nmf(d, 3, nrun = 500))
# summary(si) ## To check silhouette width
# plot(si)
# plot(si, col = c("red", "blue", "green"))
#########################################################
#########################################################
library(ConsensusClusterPlus)
dir.create("1202")
title="1202"
results = ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8,pFeature=1, title=title,
                               clusterAlg="hc",distance="pearson", seed=12345,plot="png")
############## PAC implementation ##############
maxK = 20
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i
# The optimal K
optK = Kvec[which.min(PAC)]

##########################################################
##########################################################
# Expression data
geneExp <- FPKM_UQ_Value
geneExp <- log2(geneExp + 1)

# Methylation data 
Tumor_BMIQ <- Beta_Value.impute
Tumor_BMIQ$variance <- NULL; Tumor_BMIQ$SD <- NULL

# Select subset of Pancreatic adenocarcinoma ductal

cancerExp <- geneExp[,common_sample]
cnts <- cancerExp[rowSums(cancerExp==0)< ncol(cancerExp)*0.20,] ## Remove all gene which have 25% zero's
keep <- rowSums(edgeR::cpm(cnts)>1) >= ncol(cnts)*0.20 #### Select only genes which have have CPM > 1 for >=50% samples
cnts <- cnts[keep,]
gene <- rownames(cnts)


################ MERGE METHYLATION AND EXPRESSION DATA ################################

hm450.anno <-read.table("hm450.hg38.manifest.tsv.gz",header=T,sep="\t",na.strings="",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)

rownames(hm450.anno) <- hm450.anno$probeID

#hm450.anno.unmasked <- hm450.anno[hm450.anno$MASK.general == 'FALSE',]
library(dplyr)
hm450.anno.unmasked <- hm450.anno %>% filter(MASK_general %in% "FALSE")

rm(hm450.anno)

hm450.anno.epig <-read.table("hm450.manifest.EpigeneticallySilence_Sort.tsv",header=T,sep="\t",na.strings="NA",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)

## hm450.manifest.EpigeneticallySilence_Sort.tsv is the file for CpG sites with distance from TSS. I generated this file from GDC harmonized DNA methylation file for genCode22

hm450.anno.epig <- hm450.anno.epig[which(hm450.anno.epig$Position_to_TSS!='.'),]

hm450.anno.epig.tss1500 <- hm450.anno.epig[abs(as.numeric(hm450.anno.epig$Position_to_TSS)) <= 1500,]

hm450.anno.epig.tss1500.unmasked <- hm450.anno.epig.tss1500[hm450.anno.epig.tss1500$ProbeID%in%hm450.anno.unmasked$probeID,]

hm450.tss1500 <- hm450.anno.epig.tss1500.unmasked[,c("ProbeID", "Transcript_ID", "Gene_Symbol","Position_to_TSS", "Chromosome", "Gene_Type" )]

colnames(hm450.tss1500) <- gsub("Chromosome","Chr", colnames(hm450.tss1500))

#TCGA.Harmonized <- read.table("TCGA_GDC_Harmonided_RNASeq-GFT.txt",header=T,sep="\t",na.strings="NA",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)

TCGA.Harmonized <- read.table("TCGA_GDC_Harmonided_Uniq_GTF.txt",header=T,sep="\t",na.strings="NA",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)

colnames(TCGA.Harmonized) <- c("Chr", "Reference", "Transcript", "Start", "End", "Strand", "Transcript_ID", "Gene_Symbol")

masked <- which(TCGA.Harmonized$Chr%in% c("chrX", "chrY", "chrM"))

TCGA.Har.Masked <- TCGA.Harmonized[-masked,]

TCGA.Har.Masked <- TCGA.Har.Masked[,c("Chr", "Transcript_ID", "Gene_Symbol")]

masked <- which(hm450.tss1500$Chr%in% c("chrX", "chrY", "chrM"))

hm450.tss1500.Masked <- hm450.tss1500[-masked,]

hm450.tss1500.Masked.Uniq <- unique(hm450.tss1500.Masked[,c("ProbeID", "Gene_Symbol", "Chr","Gene_Type")]) ## Unique gene and probe within +/- 1500 from TSS

# Add Ensembl Gene ID
hm450.tss1500.Masked.Uniq$gene_id <- TCGA.Harmonized$Transcript_ID[match(hm450.tss1500.Masked.Uniq$Gene_Symbol,TCGA.Harmonized$Gene_Symbol)]



# Tumor expression data
cnts_tumor <- cancerExp[gene,]
#cnts_tumor = cnts_tumor[,substr(colnames(cnts_tumor),1,15) %in%  TCGA_PAAD_select.ductal$Sample.ID]
#colnames(cnts_tumor) <- substr(colnames(cnts_tumor),1,15)

# Methylation data
#Tumor_BMIQ <- read.csv("Tumor_BMIQ.txt",header = TRUE, check.names = FALSE)
meth_data <- Tumor_BMIQ


common_patients <-Reduce(intersect,list(colnames(cnts_tumor),colnames(meth_data)))
length(common_patients)

# Take forward only tumor samples that are common in both
cnts_tumor <- cnts_tumor[,common_patients]
meth_data <- meth_data[,common_patients]
hm450.tss1500.Masked.Uniq$ENSG <- sapply(strsplit(hm450.tss1500.Masked.Uniq$gene_id, "\\."), "[[", 1)

merge1 <-cbind(hm450.tss1500.Masked.Uniq, cnts_tumor[match(hm450.tss1500.Masked.Uniq$ENSG,rownames(cnts_tumor)),])
merge1 <- na.omit(merge1)
dim(merge1)

# Merge methylation data with probe_id
merge2 <-cbind(merge1, meth_data[match(merge1$ProbeID,rownames(meth_data)),])
merge2<- na.omit(merge2)
dim(merge2)
##columns : 6 to 129 expression data 124 Cancer)
##columms :130:253 methylation data 124 cancer)
write.csv(merge2,file="exp_meth_ductal_data_11_6_2018.txt", row.names = FALSE)


#Protein Coding genes 
merge2_protein_coding <- merge2[merge2$Gene_Type=="protein_coding",]
merge2_protein_coding <- na.omit(merge2_protein_coding)
dim(merge2_protein_coding)
write.csv(merge2_protein_coding,file = "exp_meth_protein_coding_11_6_2018.txt",row.names = FALSE)

#lincRNA
merge2_lincRNA <- merge2[merge2$Gene_Type=="lincRNA",]
merge2_lincRNA <- na.omit(merge2_lincRNA)
dim(merge2_lincRNA)
write.csv(merge2_lincRNA,file = "exp_meth_lincRNA_coding_11_6_2018.txt",row.names = FALSE)

###################################################################################
###################################################################################
### Spearman correlation  
## Run this in parts : first for protein coding then for lincRNA
## Make 3 file, 1: Meth, 2: Expression, 3: Annotation
## Protein Coding
meth_file1 <- merge2_protein_coding[,93:178]
exp_file2 <- merge2_protein_coding[,7:92]
annotation_file3 <- merge2_protein_coding[,1:6]


ll <- mapply(function(x,y)cor.test(as.numeric(meth_file1[x,]),as.numeric(exp_file2[y,]), method = "spearman", alternative = "t"),
             1:nrow(meth_file1),
             1:nrow(exp_file2),
             SIMPLIFY=FALSE)

cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
spearman.p <- cbind(cor.value, p.value)
rownames(spearman.p) <- rownames(annotation_file3)
spearman.p <- data.frame(spearman.p)
spearman.p$gene <- annotation_file3$Gene_Symbol
spearman.p$gene_id <- annotation_file3$gene_id
spearman.p$Probe <- annotation_file3$ProbeID
spearman.p$chr <- annotation_file3$Chr
spearman.p <- spearman.p[!as.character(spearman.p$chr) %in% c("chrM","chrX", "chrY"),]
spearman.p$adjP <- p.adjust(spearman.p$p.value, method = c("BH"))
spearman.p1 <- spearman.p[which((spearman.p$p.value <= 0.01)&abs(spearman.p$cor.value) >=0.50),]
spearman.p1 <- spearman.p1[order(spearman.p1$gene),]
spearman.p1.neg <- spearman.p1[which(spearman.p1$cor.value < 0),]
spearman.p1.pos <- spearman.p1[which(spearman.p1$cor.value > 0),]

dim(spearman.p1.pos)
dim(spearman.p1.neg)
################################################################
write.csv(spearman.p1.neg,file="spearman_Negative_Protein_12_03_2019.txt",row.names = FALSE)
write.csv(spearman.p1.pos,file="spearman_Postive_Protein_12_03_2019.txt",row.names = FALSE)

##lincRNA
meth_file1 <- merge2_lincRNA[,93:178]
exp_file2 <- merge2_lincRNA[,7:92]
annotation_file3 <- merge2_lincRNA[,1:6]

ll <- mapply(function(x,y)cor.test(as.numeric(meth_file1[x,]),as.numeric(exp_file2[y,]), method = "spearman", alternative = "t"),
             1:nrow(meth_file1),
             1:nrow(exp_file2),
             SIMPLIFY=FALSE)

cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
spearman.p <- cbind(cor.value, p.value)
rownames(spearman.p) <- rownames(annotation_file3)
spearman.p <- data.frame(spearman.p)
spearman.p$gene <- annotation_file3$Gene_Symbol
spearman.p$gene_id <- annotation_file3$gene_id
spearman.p$Probe <- annotation_file3$ProbeID
spearman.p$chr <- annotation_file3$Chr
spearman.p <- spearman.p[!as.character(spearman.p$chr) %in% c("chrM","chrX", "chrY"),]
spearman.p$adjP <- p.adjust(spearman.p$p.value, method = c("BH"))
spearman.p1 <- spearman.p[which((spearman.p$p.value <= 0.01)&abs(spearman.p$cor.value) >=0.25),]
spearman.p1 <- spearman.p1[order(spearman.p1$gene),]
spearman.p1.neg <- spearman.p1[which(spearman.p1$cor.value < 0),]
spearman.p1.pos <- spearman.p1[which(spearman.p1$cor.value > 0),]

dim(spearman.p1.pos)

dim(spearman.p1.neg)

write.csv(spearman.p1.neg,file="spearman_Negative_lincRNA_11_6_2018.txt",row.names = FALSE)
write.csv(spearman.p1.pos,file="spearman_Postive_lincRNA_11_6_2018.txt",row.names = FALSE)



##########################################################
##########################################################
save.image("Osteosarcoma_GDC_NMF.RData")
##########################################################