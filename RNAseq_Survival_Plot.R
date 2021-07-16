setwd("D:/OneDrive - University of Nebraska Medical Center/Osteosarcoma/TARGET-OS/")
load("TARGET_OS_Survival.RData")
rm(list=setdiff(ls(), c("OS.Clinical", "OS.Clinical.TARGET", "OS.HTSeq.FPKM", "", "OS.HTSeq.FPKM",
                        "OS.HTSeq.Counts", "dataset.protein.coding", "dataset.lincRNA")))

source("plot_surv_RNAseq_function.R")


unlink("RNAseq_Surv/KM_Plot_*") 
results <- plot_surv(dir = "RNAseq_Surv", clinical_patient = PAAD.Clinical, dataGE = exp, Genelist = rownames(exp), Survresult = FALSE, Median = TRUE, p.cut = 0.05)


unlink("RNAseq_Surv/KM_Plot_*") 
results.DM <- plot_surv(dir = "RNAseq_Surv", clinical_patient = PAAD.Clinical, dataGE = exp, Genelist = rownames(results), Survresult = FALSE, Median = TRUE, p.cut = 0.05)

unlink("RNAseq_Surv/KM_Plot_*") 
results.DM.P.0.01 <- plot_surv(dir = "RNAseq_Surv", clinical_patient = PAAD.Clinical, dataGE = exp, Genelist = rownames(results.DM), Survresult = TRUE, Median = TRUE, p.cut = 0.01)


unlink("RNAseq_Surv/KM_Plot_*") 
results.DM.P.0.01 <- plot_surv(dir = "RNAseq_Surv", clinical_patient = PAAD.Clinical, dataGE = exp, Genelist = rownames(results.DM.P.0.01), Survresult = TRUE, Median = TRUE, p.cut = 0.01)

#############################################
save.image("RNASeq_survival_analysis.RData")