# R code to download GEO data for Osteosarcoma
#S. Peter Wu et al. (2017), "DNA Methylation-Based Classifier for Accurate Molecular Diagnosis of Bone Sarcomas"
#https://ascopubs.org/doi/10.1200/PO.17.00031
#JCO Precis Oncol. 2017;2017. doi: 10.1200/PO.17.00031.
library(GEOquery)
getGEOSuppFiles("GSE97529")

#Unpack the CEL files
setwd("GSE97529/")

untar("GSE97529_RAW.tar", exdir="data")
idat = list.files("data/", pattern = "idat")
sapply(paste("data", idat, sep="/"), gunzip)

txt = list.files("data/", pattern = "csv|bpm|xlsx")
sapply(paste("data", txt, sep="/"), gunzip)


gse <- getGEO("GSE97529", GSEMatrix = TRUE)
dim(pData(gse[[1]]))
head(pData(gse[[1]])[, 30:38])

write.csv(pData(gse[[1]]), "GSE97529_pData.csv")
rm(df1, idat, txt, gse)

## We don't need save anything. 
## This will download IDAT file and unzip it and save in GSE97529/data folder.