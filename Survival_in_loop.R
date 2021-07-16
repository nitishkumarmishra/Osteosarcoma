library(survival)
library(survminer)
data("lung")
#covariates <- c("age", "sex",  "ph.karno", "ph.ecog", "wt.loss")
covariates <- colnames(lung)[4:ncol(lung)]
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = lung)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("Beta", "HR (95% CI)", "Wald test", 
                                       "P-value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

univ_Cox_models <- lapply( univ_formulas, function(x){survdiff(x, data = lung)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- 1 - pchisq(x$chisq, length(x$n)-1)
                         res <- c
                         })
res.cox <- survdiff(Surv(time, status) ~ sex, data = lung)
p <- 1 - pchisq(res.cox$chisq, length(res.cox$n)-1)


df <- structure(list(times = c(724L, 1624L, 1569L, 2532L, 1271L, 2442L, 757L, 848L, 3675L, 1229L, 1582L, 1257L, 1270L, 555L, 357L, 1133L, 633L), 
                     samples = structure(c(1L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L), 
                                         .Label = c("Sample1", "Sample10", "Sample11", "Sample12", "Sample13", "Sample14", "Sample15", "Sample16", 
                                                    "Sample17", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Sample7", "Sample8", "Sample9"), 
                                         class = "factor"), 
                     vital_status = c(1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 1L), 
                     years = c(1.983561644, 4.449315068, 4.298630137, 6.936986301, 3.482191781, 6.690410959, 
                               2.073972603, 2.323287671, 10.06849315, 3.367123288, 4.334246575, 
                               3.443835616, 3.479452055, 1.520547945, 0.978082192, 3.104109589, 1.734246575), 
                     Gene1 = c(0.9, 0.8, 0.6, 1.2, 3.8, 2.3, 3.8, 0.4, 0.5, 1.2, 7.7, 2.1, 0.8, 1.8, 2.4, 3, 0.6), 
                     Gene2 = c(1.2, 3.8, 2.3, 3.8, 0.4, 0.5, 1.2, 7.7, 2.1, 0.9, 0.8, 0.6, 0.5, 1.2, 7.7, 2.1, 0.6), 
                     Gene3 = c(2.3, 3.8, 0.4, 0.5, 1.2, 7.7, 0.9, 0.8, 0.6, 0.5, 1.2, 7.7, 2.1, 0.6, 0.9, 0.8, 0.6), 
                     Gene4 = c(3.8, 0.4, 0.5, 1.2, 7.7, 2.1, 0.8, 1.8, 2.4, 3, 0.6, 0.9, 0.8, 0.6, 1.2, 3.8, 2.3), 
                     Gene5 = c(0.5, 1.2, 7.7, 0.9, 0.8, 0.6, 0.5, 1.2, 7.7, 2.1, 0.6, 0.9, 1.2, 7.7, 2.1, 0.9, 0.8)), 
                class = "data.frame", row.names = c(NA, -17L))


library(survminer)
library(survival)
# vector with the variables to run through
#genes <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5") 
genes <- names(df)[!names(df) %in% c("times",  "samples", "vital_status","years")]

for(i in 1:length(genes)){
  surv_rnaseq.cut <- surv_cutpoint(
    df,
    time = "years",
    event = "vital_status",
    variables = c(genes[i]))
  
  pdf(paste0(genes[i], "_Cuttpt.pdf"))
  print(plot(surv_rnaseq.cut, genes[i], palette = "npg"), newpage=F
  )
  dev.off()
  
  surv_rnaseq.cat <- surv_categorize(surv_rnaseq.cut)
  
  
  fit <- survfit(as.formula(paste0("Surv(years, vital_status) ~", genes[i])),
                 data = surv_rnaseq.cat)
  
  pdf(paste0(genes[i], "_Survival_high_vs_low_WithPvalue.pdf"))
  
  #plot.new() 
  print(
    ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               #linetype = "strata", # Change line type by groups
               #surv.median.line = "hv", # Specify median survival
               ggtheme = theme_survminer(), # Change ggplot2 theme
               #palette = c("#FF0027", "#060606"),
               palette = "aaas",
               xlim = c(0,10),
               break.x.by = 3,
               xlab="Time in years",
               risk.table.y.text.col = T, # colour risk table text annotations.
               risk.table.y.text = FALSE), newpage=F
  )
  
  dev.off()
}
