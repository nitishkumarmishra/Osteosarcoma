library(TCGAbiolinks)
library(SummarizedExperiment)
######################################################
## Download Clinical data from GDC #############
PRAD.Clinical <- GDCquery_clinic("TCGA-PRAD")

########### HTSeq - Counts download ##################
query.exp.hg38 <- GDCquery(project = "TCGA-PRAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - Counts")
GDCdownload(query.exp.hg38)
## Remove extra lines which is not required ###### 
PRAD.HTSeq.Counts <- GDCprepare(query = query.exp.hg38, summarizedExperiment = FALSE)
PRAD.RNASeq.Report <- PRAD.HTSeq.Counts[1:5,]
PRAD.HTSeq.Counts <- PRAD.HTSeq.Counts[-c(1:5),]
rownames(PRAD.HTSeq.Counts) <- PRAD.HTSeq.Counts$X1
PRAD.HTSeq.Counts$X1 <- NULL

########## HTSeq - FPKM download #####################
query.exp.hg38 <- GDCquery(project = "TCGA-PRAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")

GDCdownload(query.exp.hg38)
PRAD.FPKM <- GDCprepare(query = query.exp.hg38, summarizedExperiment = FALSE)
rownames(PRAD.FPKM) <- PRAD.FPKM$X1
PRAD.FPKM$X1 <-NULL


######### HTSeq - FPKM-UQ download ###################
query.exp.hg38 <- GDCquery(project = "TCGA-PRAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query.exp.hg38)
PRAD.FPKM.UQ <- GDCprepare(query = query.exp.hg38, summarizedExperiment = FALSE)
rownames(PRAD.FPKM.UQ) <- PRAD.FPKM.UQ$X1
PRAD.FPKM.UQ$X1 <- NULL

##################################################
library(dplyr)
library(DESeq2)

#####################################################
########## GTF file from GDC for ##############
#GDC.h38 GENCODE v22 GTF (used in RNA-Seq alignment and by HTSeq)
#https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
gtf <- rtracklayer::import('gencode.v22.annotation.gtf')
gtf_df=as.data.frame(gtf)
gtf_df <- gtf_df[gtf_df$type=="gene",]
gtf_df <- gtf_df[,c("gene_id", "gene_name", "gene_type")]
gtf_df <- gtf_df[,c("gene_id", "gene_name", "gene_type")]
rownames(gtf_df) <- gtf_df$gene_id

########## TCGA data structure ##########
table(substr(colnames(PRAD.HTSeq.Counts), 14, 15))
#   01    06  11 
#   498   1   52 
############# clinical file from PANCAN cancer analysis paper in Cell ##############
PANCAN.clinical <- readr::read_delim("clinical_PANCAN_patient_with_followup.tsv", delim = "\t") %>%
  filter(acronym == "PRAD") %>%
  select(bcr_patient_barcode, gender, vital_status, radiation_therapy, days_to_death, days_to_last_followup,age_at_initial_pathologic_diagnosis, race, ethnicity)  %>%
  mutate(Status=if_else(vital_status =="Alive", days_to_last_followup, days_to_death))

### This is for clinical data from cBioportal, which I downloaded. I used it for extra clinical information #
### Some time GDC have limited clinical information
cBioPortal.clinical <- readr::read_delim("prad_tcga_clinical_data.tsv", delim = "\t") %>%
  select("Patient ID" ,starts_with("Gleason"), "Sample ID") %>%
  dplyr::rename(bcr_patient_barcode="Patient ID", SampleID="Sample ID", GS_Primary="Gleason pattern primary", GS_Secondary="Gleason pattern secondary") %>%
  mutate(GS_score= GS_Primary+GS_Secondary)

#clinical <- inner_join(PANCAN.clinical, cBioPortal.clinical ## By default it will add based on common column
clinical <- inner_join(PANCAN.clinical, cBioPortal.clinical,  by = "bcr_patient_barcode")

### This is clinical file from supplementary table from TCGA PRDA paper for race information
cell.clinical <- readr::read_delim("PRAD_Cell.txt", delim = "\t") %>%
  select(PATIENT_ID ,Race) %>%
  dplyr::rename(bcr_patient_barcode=PATIENT_ID)

clinical$GS_Group <- ifelse(clinical$GS_score <= 6, "Group1", ifelse(clinical$GS_score >= 8, "Group4", ifelse(clinical$GS_Primary == 3 & clinical$GS_Secondary == 4, "Group2", "Group3")))
rownames(clinical) <- clinical$bcr_patient_barcode

write.csv(clinical, "Clinical_Data.csv")
# Few samples are duplicated, remove one of them
#TCGA-HC-7740-01 TCGA-HC-8258-01 TCGA-HC-8265-01 #These are duplicated
dataset <- PRAD.HTSeq.Counts[,grep("-01A|01B", colnames(PRAD.HTSeq.Counts))]
dataset1 <- dataset[,!grepl("TCGA-HC-7740-01B|TCGA-HC-8258-01B|TCGA-HC-8265-01B", colnames(dataset))]
colnames(dataset1) <- substr(colnames(dataset1), 1, 12)
clinical <- clinical[colnames(dataset1),]
#dataset1 <- dataset1[, as.character(PhenoData$Pheno)]
dataset1 <- dataset1[apply(dataset1,1,function(x) sum(x==0))<ncol(dataset1)*0.80,]
dataset1 <- dataset1[apply(dataset1,2,function(x) sum(x==0))<nrow(dataset1)*0.80,]

##################################################################
########################## DEG analysis ##########################
##################################################################
## Group2 vs Group1
############### Group2 vs Group1 DESeq2 analysis ########
clinical_grou1_vs2 <- clinical %>%
  filter(GS_Group=="Group1"| GS_Group=="Group2") %>%
  arrange(desc=GS_Group)
dataset1_vs2 <- dataset1[,clinical_grou1_vs2$bcr_patient_barcode]

countdata <- dataset1_vs2
condition <- factor(ifelse(clinical_grou1_vs2$GS_Group=="Group1", "control", "exp"), levels = c("control", "exp"))
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)
# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=FALSE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "gene_id"
head(resdata)
resdata <- na.omit(resdata) ## Remove genes with Pvalue = NA
resdata <- merge(gtf_df, resdata, by="gene_id")

## Write results in CSV
write.csv(clinical_grou1_vs2, "clinical_group1_vs_group2.csv")
write.csv(resdata, file="Group2_respect_Group1_diffexpr-results.csv")
resdata_0.05 <- resdata[(abs(resdata$log2FoldChange)> 1 & resdata$padj <= 0.05) ,]
write.csv(resdata_0.05, file="Group2_respect_Group1_diffexpr_log2FC_1.0_FDR_0.05.csv")
resdata_0.01 <- resdata[(abs(resdata$log2FoldChange)> 1.5 & resdata$padj <= 0.01) ,]
write.csv(resdata_0.01, file="Group2_respect_Group1_diffexpr_log2FC_1.5_FDR_0.01.csv")

