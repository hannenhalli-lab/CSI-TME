library(dplyr); library(magrittr); library(glue)
library(doParallel); library(tidyr)
library(preprocessCore); library(ggplot2)


source('../functionsForInteractionAnalysis.R')

genesets <- "development"
idh = "mut"
cohort <- "TCGA"
dims = "10"
load(glue("ICA-JADE-{dims}.signatureGenes.{tolower(cohort)}.IDH{idh}--filtered-noNorm-optimalRank.Rda"))
signatures.cellType <- signature.list[['Malignant']]

expression.file <- glue("KallistoTPMvaluesforGeneExpressionInHindbrain.txt")
sampleNames <-  data.table::fread(paste(expression.file), sep = "\t", header = F, check.names = F, skip = 0, nrow = 1) %>% as.character()
sampleNames <- sampleNames[-1]
expression.data <- read.table(expression.file, sep = "\t", header = T, check.names = F, row.names = 1, colClasses=c("character", rep("numeric",length(sampleNames))))

GeneLevelTPM <- expression.data %>% data.matrix() %>% normalize.quantiles() %>% 
  set_rownames(rownames(expression.data)) %>% 
  set_colnames(colnames(expression.data)) %>% data.frame(check.names = F)
 
 
ggTheme <- theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1), axis.text.y = element_text(size = 10),
		axis.title = element_text(size = 10, face ="bold"), panel.border = element_rect(color = "grey", fill = NA, size = 2),
		panel.grid.major = element_blank(), #panel.grid.minor.x = element_line(colour="black", size=1)
		strip.text.x = element_text(size = 10, colour = "brown")
		)


data.plot.list <- list(); target.genes <- "negative"; output.file <- glue("plots/ICA-JADE-{dims}.signatureGenes-{target.genes}.{tolower(cohort)}.IDH{idh}-expression.development.svg")
for (component in names(signatures.cellType)) {
	component.genes <- signatures.cellType[[component]][[target.genes]]
	
	if (length(component.genes) > 0) {
		data.plot <- GeneLevelTPM[rownames(GeneLevelTPM) %in% component.genes, ] %>% data.matrix %>% reshape2::melt()
		data.plot$size <- length(component.genes)
		data.plot.list[[component]] <- data.plot
	} 
}

data.plot.df <- bind_rows(data.plot.list, .id = "component")
data.plot.df$timeline <- gsub("days_Rep.+", "", data.plot.df$Var2) %>% as.numeric()
data.plot.df$stage <- ifelse(data.plot.df$timeline <= 56, "Embryonic", ifelse(data.plot.df$timeline < 280, "Fetal", "Postbirth"))
data.plot.df$expression <- log2(data.plot.df$value + 0.001)
data.plot.df <- data.plot.df %>% group_by(component, Var1, stage) %>% summarize(expression = mean(expression))
comparison <- list(c("Embryonic", "Postbirth"), c("Fetal", "Postbirth"))

box.plot <- ggplot(data = data.plot.df, aes(x = stage, y = expression)) + geom_boxplot() + facet_wrap(~component) + theme_bw() + ggTheme +
	ggpubr::stat_compare_means(label = "p.format", comparisons = comparison)


svg(file = output.file, height = 7, width = 5)
box.plot
dev.off()  

###################dimension 7 plot for negative genes##################

data.plot.list <- list(); target.genes <- "negative"; output.file <- glue("plots/ICA-JADE-{dims}.signatureGenes-{target.genes}.{tolower(cohort)}.IDH{idh}-expression.development-Dim7.svg")
for (component in names(signatures.cellType)) {
	component.genes <- signatures.cellType[[component]][[target.genes]]
	
	if (length(component.genes) > 0) {
		data.plot <- GeneLevelTPM[rownames(GeneLevelTPM) %in% component.genes, ] %>% data.matrix %>% reshape2::melt()
		data.plot$size <- length(component.genes)
		data.plot.list[[component]] <- data.plot
	} 
}

data.plot.df <- bind_rows(data.plot.list, .id = "component")
data.plot.df$timeline <- gsub("days_Rep.+", "", data.plot.df$Var2) %>% as.numeric()
data.plot.df$stage <- ifelse(data.plot.df$timeline <= 56, "Embryonic", ifelse(data.plot.df$timeline < 280, "Fetal", "Postbirth"))
data.plot.df$expression <- log2(data.plot.df$value + 0.001)
data.plot.df <- data.plot.df %>% group_by(component, Var1, stage) %>% summarize(expression = mean(expression))
data.plot.df <- data.plot.df[data.plot.df$component == "Dim.7", ]
comparison <- list(c("Embryonic", "Postbirth"), c("Fetal", "Postbirth"))

box.plot <- ggplot(data = data.plot.df, aes(x = stage, y = expression)) + geom_boxplot() + facet_wrap(~component) + theme_bw() + ggTheme +
	ggpubr::stat_compare_means(label = "p.format", comparisons = comparison)


svg(file = output.file, height = 4.5, width = 3)
box.plot
dev.off()  


