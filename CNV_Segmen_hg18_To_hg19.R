## R code for liftover from hg18 to hg19
## TARGET hg18 segmentation file to hg19
## https://www.bioconductor.org/packages/release/bioc/vignettes/genomation/inst/doc/GenomationManual.html
## https://bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html
## Further I will use it for getting SCNA data as proposed in TCGA workflow (below two URL's)
## https://bioconductor.org/packages/release/workflows/vignettes/TCGAWorkflow/inst/doc/TCGAWorkflow.html
## https://f1000research.com/articles/5-1542
###############################################################
list.files(system.file('extdata',package='genomationData'))
sampleInfo = read.table(system.file("extdata/SamplesInfo.txt", package = "genomationData"),
header = TRUE, sep = "\t")
sampleInfo[1:5, 1:5]

library(genomation)
tab.file1 = system.file("extdata/tab1.bed", package = "genomation")
readGeneric(tab.file1, header = TRUE)

###### I just make changes in tab2.bed which similar to my CNV segmentation file
#start	end	chr	count	sum	avg	Target
#9439473	9437272	chr21	285	1426	25.9	TARGET-40-0A4HLD-01A-01D
#9484663	9483485	chr21	165	818	28	TARGET-40-0A4HLD-01A-01D
#9648116	9647866	chr21	18	168	14.4	TARGET-40-0A4HLD-01A-01D
#9709231	9708935	chr21	31	218	20.9	TARGET-40-0A4HLD-01A-01D
#9826296	9825442	chr21	120	568	28.1	TARGET-40-0A4HLD-01A-01D
#9909218	9909011	chr21	20	143	19.3	TARGET-40-0A4HLD-01A-01D

###############################################################
library(rtracklayer)

tab.file2 = system.file("extdata/tab2.bed", package = "genomation")
file <- readGeneric(tab.file2, chr = 3, start = 2, end = 1, keep.all.metadata = TRUE, header = TRUE)
genome(file) = "hg18"

#path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
path = system.file(package="liftOver", "extdata", "hg18ToHg19.over.chain")
#path = system.file(package="liftOver", "data", "hg18ToHg19.over.chain")
ch = import.chain(path)

seqlevelsStyle(file) = "UCSC"  # necessary
cur19 = liftOver(file, ch)
class(cur19)

cur19 = unlist(cur19)
genome(cur19) = "hg19"
#cur19 = new("gwaswloc", cur19)
cur19
file
#############################################################
#############################################################
setwd(dir = "D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma")

#head(read.csv("TARGET_OS_L3_Segmentation.txt", header = TRUE, sep = "\t", check.names = FALSE))
segmentation <- readGeneric(file = "TARGET_OS_L3_Segmentation.txt", chr = 2, start = 3, end = 4, keep.all.metadata = TRUE, header = TRUE)
genome(segmentation) = "hg18"

# put SampleName at last
# use plyranges which will handle Granges as data.frame, I can use dplyr wrangling command on Granges data.
library(plyranges)
segmentation <- segmentation %>%
  select(-SampleName, everything())


library(rtracklayer)
#path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
path = system.file(package="liftOver", "extdata","hg18ToHg19.over.chain")
ch = import.chain(path)

seqlevelsStyle(segmentation) = "UCSC"  # necessary
segmentation19 = liftOver(segmentation, ch)
class(segmentation19)

segmentation19 = unlist(segmentation19)
genome(segmentation19) = "hg19"
#cur19 = new("gwaswloc", cur19)
segmentation19
segmentation


save.image("CNV_Segmen_hg18_To_hg19.RData")

#########################################################
#########################################################
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-==--=--==---=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-==--=--==--
# Retrieve probes meta file from broad institute website for hg19
# For hg38 analysis please take a look on:
# https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files
# File: SNP6 GRCh38 Liftover Probeset File for Copy Number Variation Analysis
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-==--=--==---=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-==--=--==--
library(TCGAbiolinks)
library(downloader)
library(readr)
library(gaia)
library(tidyverse)
# Select common CN technology available for GBM and LGG


df1 <- data.frame(iranges = segmentation19)
colnames(df1) <- gsub("iranges.", "", colnames(df1))
#df3 <- data.frame(iranges = segmentation)
df2 <- df1 %>%
  select(SampleName, seqnames, start,  end, Num.Probes, Segment.Mean) %>%
  #select(Sample = SampleName, Chromosome = seqnames, Start = start,  End = end, Num_Probes = Num.Probes, Segment_Mean = Segment.Mean) %>%
  mutate(seqnames = gsub("chr", "", seqnames)) %>%
  rename(Sample = SampleName,  Chromosome = seqnames, Start = start,  End = end, Num_Probes = Num.Probes, Segment_Mean = Segment.Mean) %>%
  as_tibble()

##############################################################
##############################################################
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-==--=--==---=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-==--=--==--
gdac.root <- "ftp://ftp.broadinstitute.org/pub/GISTIC2.0/hg19_support/"
file <- paste0(gdac.root, "genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1.txt")
if(!file.exists(basename(file))) downloader::download(file, basename(file))
markersMatrix <-  readr::read_tsv(basename(file), col_names = FALSE, col_types = "ccn", progress = FALSE)
#save(markersMatrix, file = "markersMatrix.rda", compress = "xz")

save(markersMatrix, file = "markersMatrix.rda", compress = "xz")

#########################################################################
#########################################################################
###############################
## CNV  data  pre-processing ##
###############################
cancer <- "TARGET-OS"
message(paste0("Starting ", cancer))
# get objects created above
#data(GBMnocnvhg19)
#data(markersMatrix)
cnvMatrix <- df2

# Add label (0 for loss, 1 for gain)
cnvMatrix <- cbind(cnvMatrix,Label=NA)
cnvMatrix[cnvMatrix[,"Segment_Mean"] < -0.3,"Label"] <- 0
cnvMatrix[cnvMatrix[,"Segment_Mean"] > 0.3,"Label"] <- 1
cnvMatrix <- cnvMatrix[!is.na(cnvMatrix$Label),]

# Remove "Segment_Mean" and change col.names
cnvMatrix <- cnvMatrix[,-6]
colnames(cnvMatrix) <- c("Sample.Name", "Chromosome", "Start", "End", "Num.of.Markers", "Aberration")

# Substitute Chromosomes "X" and "Y" with "23" and "24"
cnvMatrix[cnvMatrix$Chromosome == "X","Chromosome"] <- 23
cnvMatrix[cnvMatrix$Chromosome == "Y","Chromosome"] <- 24
cnvMatrix$Chromosome <- as.integer(cnvMatrix$Chromosome)

# Recurrent CNV identification with GAIA
colnames(markersMatrix) <- c("Probe.Name", "Chromosome", "Start")
unique(markersMatrix$Chromosome)
######################################################################
######################################################################
markersMatrix[markersMatrix$Chromosome == "X","Chromosome"] <- "23"
markersMatrix[markersMatrix$Chromosome == "Y","Chromosome"] <- "24"
markersMatrix$Chromosome <- as.integer(markersMatrix$Chromosome)
markerID <- paste(markersMatrix$Chromosome,markersMatrix$Start, sep = ":")
# Removed duplicates
markersMatrix <- markersMatrix[!duplicated(markerID),]
# Filter markersMatrix for common CNV
markerID <- paste(markersMatrix$Chromosome,markersMatrix$Start, sep = ":")

file <- "ftp://ftp.broadinstitute.org/pub/GISTIC2.0/hg19_support/CNV.hg19.bypos.111213.txt"
if(!file.exists(basename(file))) downloader::download(file, basename(file))
commonCNV <- readr::read_tsv(basename(file), progress = FALSE)
#data(CNV.hg19.bypos.111213)
commonID <- paste(commonCNV$Chromosome,commonCNV$Start, sep = ":")
markersMatrix_fil <- markersMatrix[!markerID %in% commonID,]


# set.seed(200)
# markers_obj <- load_markers(as.data.frame(markersMatrix_fil))
# nbsamples <- length(unique(cnvMatrix$Sample))
# cnv_obj <- load_cnv(cnvMatrix, markers_obj, nbsamples)

library(gaia)
set.seed(200)
markers_obj <- load_markers(
  marker_matrix = markersMatrix_fil %>% as.data.frame()
)

nbsamples <- length(unique(cnvMatrix$Sample.Name))
cnv_obj <- load_cnv(
  segmentation_matrix = cnvMatrix  %>% as.data.frame(), 
  markers_list = markers_obj, 
  num_of_samples = nbsamples
)


suppressWarnings({
  results <- runGAIA(cnv_obj,
                     markers_obj,
                     output_file_name = paste0("GAIA_",cancer,"_flt.txt"),
                     aberrations = -1,  # -1 to all aberrations
                     #chromosomes = 9, # -1 to all chromosomes
                     chromosomes = -1, # -1 to all chromosomes
                     approximation = TRUE, # Set to TRUE to speed up the time requirements
                     num_iterations = 5000, # Reduced to speed up the time requirements
                     threshold = 0.25)
})
# Set q-value threshold
# Use a smalled value for your analysis. We set this as high values
# due to the small number of samples which did not reproduced
# results with smaller q-values
threshold <- 0.3

# Plot the results
RecCNV <- t(apply(results,1,as.numeric))
colnames(RecCNV) <- colnames(results)
RecCNV <- cbind(RecCNV, score = 0)
minval <- format(min(RecCNV[RecCNV[,"q-value"] != 0,"q-value"]), scientific = FALSE)
minval <- substring(minval,1, nchar(minval) - 1)
RecCNV[RecCNV[,"q-value"] == 0,"q-value"] <- as.numeric(minval)
RecCNV[,"score"] <- sapply(RecCNV[,"q-value"],function(x) -log10(as.numeric(x)))
RecCNV[RecCNV[,"q-value"] == as.numeric(minval),]
gaiaCNVplot(RecCNV,threshold)
save(results, RecCNV, threshold, file = paste0(cancer,"_CNV_results.rda"))
################################################################################
#################### Gene annotation of recurrent CSV ##########################
################################################################################
library(GenomicRanges)
# Get gene information from GENCODE using biomart
genes <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg19") 
genes <- genes[genes$external_gene_name != "" & genes$chromosome_name %in% c(1:22,"X","Y"),]
genes[genes$chromosome_name == "X", "chromosome_name"] <- 23
genes[genes$chromosome_name == "Y", "chromosome_name"] <- 24
genes$chromosome_name <- sapply(genes$chromosome_name,as.integer)
genes <- genes[order(genes$start_position),]
genes <- genes[order(genes$chromosome_name),]
genes <- genes[,c("external_gene_name", "chromosome_name", "start_position","end_position")]
colnames(genes) <- c("GeneSymbol","Chr","Start","End")
genes_GR <- makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
save(genes_GR,genes,file = "genes_GR.rda", compress = "xz")

##############################
## Recurrent CNV annotation ## 
##############################
# Get gene information from GENCODE using biomart
#data(genes_GR) # downloaded in the previous step (available in TCGAWorkflowData)

load(paste0(cancer,"_CNV_results.rda"))
sCNV <- RecCNV[RecCNV[,"q-value"] <= threshold,c(1:4,6)]
sCNV <- sCNV[order(sCNV[,3]),]
sCNV <- sCNV[order(sCNV[,1]),]
colnames(sCNV) <- c("Chr","Aberration","Start","End","q-value")
sCNV_GR <- makeGRangesFromDataFrame(sCNV,keep.extra.columns = TRUE)

hits <- findOverlaps(genes_GR, sCNV_GR, type = "within")
sCNV_ann <- cbind(sCNV[subjectHits(hits),],genes[queryHits(hits),])
AberrantRegion <- paste0(sCNV_ann[,1],":",sCNV_ann[,3],"-",sCNV_ann[,4])tibble::tribble(
  ~`ECD.RNA-Seq`, ~`2_RNAseq`,  ~ECD, ~`0.272832098`, ~`-0.303806097`, ~`-0.166770617`, ~`-0.198378325`, ~`-0.122169005`, ~`-0.323743922`, ~`-0.433584441`, ~`0.104662408`, ~`0.046847872`, ~`-0.198584014`, ~`-0.017782583`, ~`-0.056533299`, ~`-0.20269372`, ~`-0.495141359`, ~`0.048527427`, ~`0.320579947`, ~`-0.659106058`, ~`-0.068965075`, ~`0.111484113`, ~`-0.122598668`, ~`0.447398006`, ~`-0.060088513`, ~`0.156971537`, ~`-0.143731949`, ~`-0.386653991`, ~`0.607736088`, ~`0.106627293`, ~`0.077753661`, ~`0.215981613`, ~`-0.154161966`, ~`-0.852205909`, ~`-0.041737235`, ~`0.162568279`, ~`0.280327122`, ~`0.125512955`, ~`0.612406386`, ~`-0.008818325`, ~`0.053109956`, ~`-0.039467691`, ~`-0.099021433`, ~`-0.128406097`, ~`0.013201161`, ~`0.145689288`, ~`0.423139419`, ~`0.424820293`, ~`0.264640205`, ~`0.032942874`, ~`0.313880518`, ~`-0.126213149`, ~`0.342526799`, ~`0.208119457`, ~`-0.04899088`, ~`0.159652931`, ~`-0.41176893`, ~`-0.020204293`, ~`-0.362874265`, ~`-0.265587329`, ~`0.130366711`, ~`0.36895247`, ~`-0.147531194`, ~`0.438152564`, ~`0.278989542`, ~`-0.262663754`, ~`-0.197847542`, ~`0.115271058`, ~`0.513072975`, ~`-0.078815639`, ~`-0.161448274`, ~`0.101806413`, ~`0.170593912`, ~`-0.014348513`, ~`0.17902147`, ~`0.14687162`, ~`-0.267018551`, ~`0.023795994`, ~`0.135394095`, ~`0.354299692`,
  "ECD MS Protein", "3_Protein", "ECD",   -0.043952866,     0.639717194,    -0.598166939,     0.314205905,       0.8558184,     -0.51044051,    -0.519710523,    1.195314303,    1.051879696,     0.418781311,    -0.450903024,    -0.395103448,   -0.252940516,    -0.907368024,   -0.745750313,    0.761815488,     0.093982707,    -0.781776657,    0.514835901,     0.078212036,    1.720602277,    -0.369478609,    -0.58018001,    -0.529029782,    -0.127706005,   -0.587412721,   -0.543410479,    1.040919049,    0.322898189,     0.285540304,     0.128116222,    -0.990971004,    0.327146495,    0.751728238,    0.291373908,   -1.267242053,    -0.308847158,   -0.179494109,    -0.024277472,     -0.72804512,    -3.070664084,    1.227991859,   -2.134612563,    0.721736931,    0.445093185,    -0.37821446,    0.023673384,    1.176838596,     0.646915082,   -0.340301326,    0.123491893,    0.308157219,    1.197600592,   -0.419892419,    -0.088969922,     0.478649142,     0.614950022,   -0.247504497,   0.708685898,      0.23854713,    0.177628256,    1.273202172,    -0.701935845,     1.372054954,    0.421987498,    0.465010756,    -1.158220944,     0.058520773,   -0.637656256,    0.063414224,     -0.30126394,   0.953675561,  -0.837818361,     1.034503258,    0.414413708,   -0.272911979,    1.477777003
)

GeneRegion <- paste0(sCNV_ann[,7],":",sCNV_ann[,8],"-",sCNV_ann[,9])
AmpDel_genes <- cbind(sCNV_ann[,c(6,2,5)],AberrantRegion,GeneRegion)
AmpDel_genes[AmpDel_genes[,2] == 0,2] <- "Del"
AmpDel_genes[AmpDel_genes[,2] == 1,2] <- "Amp"
rownames(AmpDel_genes) <- NULL