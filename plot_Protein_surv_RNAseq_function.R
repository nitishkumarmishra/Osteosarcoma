library(survival)
library(survminer)
plot_surv.Protein <- function (dir, clinical_patient, dataGE, Genelist, Survresult = FALSE, Median = TRUE, Mean =TRUE, Sdev=0, p.cut = 0.05) 
{
  if (!dir.exists(dir)) dir.create(dir)
  if (!is.null(Genelist)) {
    Genelist <- intersect(rownames(dataGE), Genelist)
    group2 <- colnames(dataGE[,grep("-01", colnames(dataGE))])
    dataCancer <- dataGE[Genelist, group2, drop = FALSE]
    colnames(dataCancer) <- substr(colnames(dataCancer), 1, 16)
  }
  cfu <- OS.Clinical.TARGET[ OS.Clinical.TARGET[, "TARGET USI"] %in% substr(colnames(OS.HTSeq.FPKM), 1, 16), ]
  cfu <- subset(cfu, select=c("TARGET USI", "Overall Survival Time in Days", "Event Free Survival Time in Days", "Vital Status"))
  colnames(cfu) <- c("bcr_patient_barcode", "days_to_death", "days_to_last_follow_up", "vital_status")
  #cfu <- as.data.frame(subset(cfu, select = c("bcr_patient_barcode", "days_to_death", "days_to_last_follow_up", "vital_status")))
  cfu$time <- cfu$days_to_death
  cfu$days_to_death <- as.numeric(as.character(cfu$days_to_death))
  cfu$days_to_last_follow_up <- as.numeric(as.character(cfu$days_to_last_follow_up))
  rownames(cfu) <- cfu[, "bcr_patient_barcode"]
  
  cfu <- cfu[!is.na(cfu$days_to_death),]
  cfu <- cfu[which(cfu$days_to_death>=10),]
  
  cfu_complete <- cfu
  
  dataCancer <- dataCancer[,rownames(cfu_complete)]
  ngenes <- nrow(dataCancer)
  
  tabSurv_Matrix <- matrix(0, nrow(as.matrix(rownames(dataCancer))), 13)
  colnames(tabSurv_Matrix) <- c("mRNA", "pvalue", "Cancer Deaths", "Cancer Deaths with Top", "Cancer Deaths with Down", "Surv Mean Top Months", "Surv Mean Down Months", "HR", "HR_low95", "HR_up95", "P_value", "beta", "Wald.test")
  tabSurv_Matrix <- as.data.frame(tabSurv_Matrix)
  
  for (i in 1:ngenes) {
    cat(paste0((ngenes - i), "."))
    mRNAselected <- rownames(dataCancer)[i]
    mRNAselected_values <- dataCancer[rownames(dataCancer) == mRNAselected, ]
    
    if (all(mRNAselected_values == 0)) 
      next
    tabSurv_Matrix[i, "mRNA"] <- mRNAselected
    mRNAselected_values_ordered <- sort(mRNAselected_values, decreasing = TRUE)
    
    if(Median==TRUE){
      sd <- sd(as.numeric(mRNAselected_values_ordered))
      mRNAselected_values_ordered_top <- median(as.numeric(mRNAselected_values_ordered))+ (sd*Sdev)
      mRNAselected_values_ordered_down <- median(as.numeric(mRNAselected_values_ordered))- (sd*Sdev)
    }
    
    if(Mean==TRUE){
      sd <- sd(as.numeric(mRNAselected_values_ordered))
      mRNAselected_values_ordered_top <- mean(as.numeric(mRNAselected_values_ordered))+ (sd*Sdev)
      mRNAselected_values_ordered_down <- mean(as.numeric(mRNAselected_values_ordered))- (sd*Sdev)
    }
    
    mRNAselected_values_newvector <- mRNAselected_values
    
    if (!is.na(mRNAselected_values_ordered_top)) {
      numberOfSamples <- length(mRNAselected_values_ordered)
      lastelementTOP <- max(which(mRNAselected_values_ordered >= mRNAselected_values_ordered_top))
      firstelementDOWN <- min(which(mRNAselected_values_ordered <= mRNAselected_values_ordered_down))
      samples_top_mRNA_selected <- names(mRNAselected_values_ordered[1:lastelementTOP])
      samples_down_mRNA_selected <- names(mRNAselected_values_ordered[firstelementDOWN:numberOfSamples])
      samples_UNCHANGED_mRNA_selected <- names(mRNAselected_values_newvector[which((mRNAselected_values_newvector) > 
                                                                                     mRNAselected_values_ordered_down & mRNAselected_values_newvector < 
                                                                                     mRNAselected_values_ordered_top)])
      cfu_onlyTOP <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_top_mRNA_selected, ]
      cfu_onlyDOWN <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_down_mRNA_selected, ]
      cfu_onlyUNCHANGED <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_UNCHANGED_mRNA_selected, ]
      cfu_ordered <- NULL
      cfu_ordered <- rbind(cfu_onlyTOP, cfu_onlyDOWN)
      cfu <- cfu_ordered
      ttime <- cfu$time
      #sum(status <- ttime > 0)
      deads_complete <- sum(cfu$vital_status=="Dead")
      deads_top <- sum(cfu_onlyTOP$vital_status=="Dead")
      
      if (dim(cfu_onlyDOWN)[1] >= 1) {deads_down <- sum(cfu_onlyDOWN$vital_status=="Dead")}
      else {deads_down <- 0 }
      tabSurv_Matrix[i, "Cancer Deaths"] <- deads_complete
      tabSurv_Matrix[i, "Cancer Deaths with Top"] <- deads_top
      tabSurv_Matrix[i, "Cancer Deaths with Down"] <- deads_down
      #tabSurv_Matrix[i, "Mean Normal"] <- mean(as.numeric(mRNAselected_values_normal))
      dataCancer_onlyTop_sample <- dataCancer[, samples_top_mRNA_selected, drop = FALSE]
      dataCancer_onlyTop_sample_mRNAselected <- dataCancer_onlyTop_sample[rownames(dataCancer_onlyTop_sample) == mRNAselected, ]
      dataCancer_onlyDown_sample <- dataCancer[, samples_down_mRNA_selected, drop = FALSE]
      dataCancer_onlyDown_sample_mRNAselected <- dataCancer_onlyDown_sample[rownames(dataCancer_onlyDown_sample) == mRNAselected, ]
      #tabSurv_Matrix[i, "Mean Tumor Top"] <- mean(as.numeric(dataCancer_onlyTop_sample_mRNAselected))
      #tabSurv_Matrix[i, "Mean Tumor Down"] <- mean(as.numeric(dataCancer_onlyDown_sample_mRNAselected))
      ttime <- cfu$time
      ttime <- ttime*0.032854884083862 ## conver days in months
      tabSurv_Matrix[i, "Surv Mean Top Months"] <- mean(ttime[1:nrow(cfu_onlyTOP)])
      tabSurv_Matrix[i, "Surv Mean Down Months"] <- mean(tail(ttime, nrow(cfu_onlyDOWN)))
      ttime1 <- ttime
      status <- ifelse(cfu$vital_status=="Dead", 1, 0)
      ttime <- Surv(ttime, status)
      rownames(ttime) <- rownames(cfu)
      legendHigh <- paste0("High", " (",nrow(cfu_onlyTOP),  " ,", deads_top,")")
      legendLow <- paste0("Low", " (",nrow(cfu_onlyDOWN), " ,", deads_down,")")
      #print(paste0("Now running analysis for the mRNA :: ",mRNAselected))
      tabSurv_pvalue <- tryCatch({
        tabSurv <- survdiff(ttime ~ c(rep("top", nrow(cfu_onlyTOP)), rep("down", nrow(cfu_onlyDOWN)))) 
        tabSurv_chis <- unlist(tabSurv)$chisq
        pvalue_survdiff <- as.numeric(1 - pchisq(abs(tabSurv$chisq), df = 1))
        fac <- c(rep("top", nrow(cfu_onlyTOP)), rep("down", nrow(cfu_onlyDOWN)))
        fac <- ifelse(fac=="top", 2, 1)## so low will be reference
        data1 <- as.data.frame(cbind(ttime1, status, fac))
        tabSurv <- survminer::surv_fit(Surv(ttime1, status) ~ fac, data = data1)
        tabSurv <- survminer::surv_pvalue(tabSurv, data = data1)
        pvalue <- NULL
        pvalue <-  tabSurv[2]
        pvalue <- round(pvalue, 3) 
        pvalue1 <- pvalue
        pvalue <- ifelse(pvalue==0, paste0("< 0.001"), pvalue)
        print(paste0("Now running for the mRNA :: ",mRNAselected, "  Survdiff P-value = ", pvalue_survdiff, "  survminer P-value = ", pvalue))
        #p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
        #pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
        x <- coxph(ttime ~ c(rep("top", nrow(cfu_onlyTOP)), rep("down", nrow(cfu_onlyDOWN))))
        x <- summary(x)
        Pvalue<-signif(x$waldtest["pvalue"], digits=2)
        wald.test<-signif(x$waldtest["test"], digits=2)
        beta<-signif(x$coefficients[1], digits=2);#coeficient beta
        HR <-signif(x$coefficients[2], digits=2);#exp(beta)
        HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
        HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
      }, error = function(e) {
        return(Inf)
      })
      tabSurv_Matrix[i, "pvalue"] <- pvalue
      tabSurv_Matrix[i, "beta"] <- beta
      tabSurv_Matrix[i, "HR"] <- HR
      tabSurv_Matrix[i, "Wald.test"] <- wald.test
      tabSurv_Matrix[i, "P_value"] <- Pvalue
      tabSurv_Matrix[i, "HR_low95"] <- HR.confint.lower
      tabSurv_Matrix[i, "HR_up95"] <- HR.confint.upper
      if (Survresult == TRUE & pvalue1 <= p.cut) {
        titlePlot <- paste("Kaplan-Meier survival analysis, pvalue = ", pvalue, "\nCox regression, ", "pvalue = ", Pvalue, ", HR = ", paste0(HR, "(CI ",HR.confint.lower, " - ", HR.confint.upper, ")"))
        par( mar=c(4.3, 4.3, 4, 1))
        # To make High as first use A for high and B for low. This way low will became reference in plot. 
        plot(survfit(ttime ~ c(rep("B", nrow(cfu_onlyTOP)), rep("A", nrow(cfu_onlyDOWN)))), mark.time=TRUE, col = c("blue4", "red4"),  main =  titlePlot, xlab = "Overall Survival Months", ylab = "Survival Probability", lwd = 2.5, cex.main=1.4, cex.lab=1.4, cex.axis=1.4)
        text(60.5, 0.98, paste0("Gene (",mRNAselected,")"), cex = 1.5, col = "blue4")
        legend("topright", bty="n", x.intersp=0.1, bg='lightblue', inset = 0, y.intersp=0.75, lty=1, lwd=1, cex=1.4, text.font=1, legend = c(legendHigh, legendLow), col = c("red4","blue4"), text.col = c("red4","blue4"), pch = 15)
        #f.horlegend("topleft", legend = c(legendHigh, legendLow), lty=1, lwd=1, text.cex = 1.4, bbord =  c("red4","blue4"), text.col = c("red4","blue4"), pch = 15)
        #print(tabSurv)
        #file_name = paste0("Survival_Plot_DM_Hypo/KM_Plot_", mRNAselected, ".pdf")
        file_name = paste0(dir, "/", mRNAselected, ".pdf")
        dev.print(pdf, file_name, width = 9, height = 9)
      }
    }
  }
  tabSurv_Matrix[tabSurv_Matrix == "-Inf"] <- 0
  tabSurvKM <- tabSurv_Matrix
  tabSurvKM <- tabSurvKM[tabSurvKM$mRNA != 0, ]
  tabSurvKM <- tabSurvKM[tabSurvKM$pvalue <= p.cut, ]
  tabSurvKM <- tabSurvKM[!duplicated(tabSurvKM$mRNA), ]
  rownames(tabSurvKM) <- tabSurvKM$mRNA
  tabSurvKM <- tabSurvKM[, -1]
  tabSurvKM <- tabSurvKM[order(tabSurvKM$pvalue, decreasing = FALSE), ]
  return(tabSurvKM)
}
#####################