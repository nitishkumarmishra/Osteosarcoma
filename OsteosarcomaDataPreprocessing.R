library(sesame)
library(ChAMP)
#GEOquery::getGEOSuppFiles("GSE72872") # R command to download IDAT files. IDAT are in supplementary data.
setwd("D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/")
# searchIDATprefixes will work in previous directory. We put all IDAT files in one directory and provide it's name.
IDATprefixes <- searchIDATprefixes("All IDAT data")
betas <- openSesame(IDATprefixes)

##########################################
# DNA methylation processing with sesame #
##########################################
## Remove CpGs which have NA for more than 20% samples
## Then do imputation by using imputeKNN
Meth <- betas
numNAs <- rowSums(is.na(Meth))
Meth <- Meth[!(numNAs > ncol(Meth)*0.2),]
Meth <- impute::impute.knn(Meth, k=15, rng.seed = 123)
Meth <- Meth$data
dir = "CHAMP_Normalization"
if (!dir.exists(dir)) unlink(dir, recursive = TRUE)
myNorm <- champ.norm(beta=Meth, method="BMIQ",
                     plotBMIQ=FALSE, arraytype="450K", cores=3)

#########################################
###### Change the name of myNorm ########
#########################################
MANIFEST <- read.csv("MANIFEST_ALL_FOR_R.txt", header = TRUE, sep = "\t")
MANIFEST.select <- MANIFEST[, c("Extract.Name","Array.Data.File")]
MANIFEST.select$ID <- gsub("_Grn.idat|_Red.idat", "", MANIFEST.select$Array.Data.File)
MANIFEST.select$Array.Data.File <- NULL
MANIFEST.select <- unique(MANIFEST.select)

GSE97529_pData <- read.csv("GSE97529/GSE97529_pData.csv", header = TRUE, sep = ",", row.names = 1)
GSE97529_pData.select <- GSE97529_pData[,c("geo_accession", "source_name_ch1", "supplementary_file")]
GSE97529_pData.select$source_name_ch1 <- gsub(" ", "-", GSE97529_pData.select$source_name_ch1)
GEO_ID <- stringr::str_split_fixed(GSE97529_pData.select$supplementary_file, pattern = "[/]", 9)[,9]
GEO_ID <- stringr::str_sub(GEO_ID, end = -13)
GSE97529_pData.select$ID <- GEO_ID
GSE97529_pData.select$supplementary_file <- NULL

rownames(MANIFEST.select) <- MANIFEST.select$Extract.Name
rownames(GSE97529_pData.select) <- GSE97529_pData.select$source_name_ch1

tt <- MANIFEST.select ; tt$Extract.Name <- NULL; tt1 <- GSE97529_pData.select
tt1$geo_accession <- NULL; tt1$source_name_ch1 <- NULL; List <- rbind(tt, tt1)
List$Name <- rownames(List); rownames(List) <- List$ID
List <- List[colnames(myNorm),]; List$ID <- NULL
colnames(myNorm) <- List$Name

rm(tt, tt1, dir, GEO_ID, numNAs)

#########################################
save.image("OsteosarcomaData.RData")
#########################################
saveRDS(myNorm, "OsteoMeth.rds")
