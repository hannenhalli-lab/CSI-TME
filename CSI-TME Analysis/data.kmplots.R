library(glue)
library(survival)
library(survminer)
source('../functionsForInteractionAnalysis.R')


cancerType <- c("GBM", "LGG")
clinical <- read.table("AllAvialableClinicalDataMod.txt", sep = "\t", header = T)
clinical <- clinical[clinical$type %in% cancerType, ]


create_km_plot <- function(clinical.data, cancer_type, groups.data, variable.name, make_groups) {

  clinical.cancer <- clinical.data[clinical.data$type %in% cancer_type, ]
  clinical_matched <- merge(clinical.cancer, groups.data, by = "samples")
  clinical_matched$group <- clinical_matched[ ,variable.name]
  
  # survival analysis for KM plot
  fit <- survfit(Surv(time, status) ~ group, data = clinical_matched)
  
  levels.names <- fit$strata %>% names %>% gsub("group=", "", .)
  
  p <- ggsurvplot(fit, data = clinical_matched, pval = TRUE, 
                  title = paste(variable.name),
                  xlab = "Time (days)", ylab = "Survival Probability",
                  legend.title = "Marker_score",
                  legend.labs = levels.names)
  
  return(p)
  
}

sampleType = "mut"
cohort <- "TCGA"; cohort.names <- "tcga"
method <- "JADE"; dims = "10"

load(glue("ICA-{method}-{dims}.cellTypes.{tolower(cohort)}.IDH{sampleType}-filtered-noNorm-optimalRank.Rda"))

binMapList <- componentBinningFunction(pcaList = ica.res.list, factorization = "ICA", sample.names = NULL)

#cellType1 <- "Malignant"; cellType2 <- "Endothelial"
#dim1 <- "Dim.2"; dim2 <- "Dim.3";
#direction1 <- "positive"; direction2 <- "negative"
#status.lr <- "deactivation"

#int.bin <- "Bin.9"

#name.ic1 <- "tip_like_endothelial"
#name.ic2 <- "EMT.Stress_Malignant"

cellType1 <- "Malignant"; cellType2 <- "Tcell"
dim1 <- "Dim.3"; dim2 <- "Dim.2";
direction1 <- "negative"; direction2 <- "negative"
int.bin <- "Bin.1"

status.lr <- "activation"


name.ic1 <- "IC3.Malignant"
name.ic2 <- "Adhesion_Tcells"


out_ic1 <- glue("plots/survival.{name.ic1}.svg")
out_ic2 <- glue("plots/survival.{name.ic2}.svg")
out_int <- glue("plots/survival.interaction.{name.ic1}-{name.ic2}-{int.bin}.svg")



if (direction1 == "positive") {
	ic1.level <- 2
}
if (direction1 == "negative") {
	ic1.level <- 0
}

if (direction2 == "positive") {
	ic2.level <- 2
}
if (direction2 == "negative") {
	ic2.level <- 0
}


comp1.bins <- binMapList[[cellType1]][ ,dim1]
comp2.bins <- binMapList[[cellType2]][ ,dim2]

ic1.levels <- ifelse(comp1.bins == ic1.level, "high", "low")
ic2.levels <- ifelse(comp2.bins == ic2.level, "high", "low")

ic.levels.df <- data.frame(ic1.levels, ic2.levels) %>% set_names(c(name.ic1, name.ic2))
ic.levels.df$samples <- gsub("\\.\\w+\\.\\w+\\.\\w+\\.\\w+$", "", rownames(ic.levels.df))
ic.levels.df$samples <- gsub("\\.", "-", ic.levels.df$samples)


#binValues <- data.frame(barcodes = names(comp1.bins), bin1 = comp1.bins, bin2 = comp2.bins) %>% set_rownames(NULL)
#binValues$samples <- gsub("\\.\\w+\\.\\w+\\.\\w+\\.\\w+$", "", binValues$barcodes)
#binValues$samples <- gsub("\\.", "-", binValues$samples)

if (int.bin == "Bin.1") {
	interaction <- ifelse(binValues$bin1 == 0 & binValues$bin2 == 0, 1, 0)
}
if (int.bin == "Bin.3") {
	interaction <- ifelse(binValues$bin1 == 2 & binValues$bin2 == 0, 1, 0)
}
if (int.bin == "Bin.9") {
	interaction <- ifelse(binValues$bin1 == 2 & binValues$bin2 == 2, 1, 0)
}
ic.levels.df$interaction <- interaction
#binValues$interaction <- interaction
#binValues <- binValues[ ,c('samples', 'bin1', 'bin2', 'interaction')]




p1 <- create_km_plot(clinical.data = clinical, cancer_type = c("LGG", "GBM"), groups.data = ic.levels.df, variable.name = name.ic1, make_groups = F)

p2 <- create_km_plot(clinical.data = clinical, cancer_type = c("LGG", "GBM"), groups.data = ic.levels.df, variable.name = name.ic2, make_groups = F)

p3 <- create_km_plot(clinical.data = clinical, cancer_type = c("LGG", "GBM"), groups.data = ic.levels.df, variable.name = "interaction", make_groups = F)

ggsave(out_ic1, p1$plot, height = 4, width = 4)
ggsave(out_ic2, p2$plot, height = 4, width = 4)
ggsave(out_int, p3$plot, height = 4, width = 4)


