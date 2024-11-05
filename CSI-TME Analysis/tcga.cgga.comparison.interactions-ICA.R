#library(clusterProfiler)
library(dplyr); library(magrittr);
library(ggplot2); library(glue);
library(purrr); library(survival);
library(ggpubr)

source('../functionsForInteractionAnalysis.R')


####do a fisher test####
fisherTestFunction <- function(list1, list2, universe) {
	
	if (is.numeric(universe)) {
		universe = universe
	} else {
		list1 <- list1[list1 %in% universe]
		list2 <- list2[list2 %in% universe]
		universe <- length(universe)
	}
	
	c1 <- sum(list1 %in% list2)
	c2 <- length(list1) - c1
	c3 <- length(list2) - c1
	c4 <- universe - (c3 + length(list1))
	
	fisherMatrix <- matrix(c(c1,c2,c3,c4), nrow = 2)
	fisherTest <- fisher.test(fisherMatrix)
	enrichment <- fisherTest$estimate
	pvalue <- fisherTest$p.value
	CI <- fisherTest$conf.int
	
	return(c(L1 = length(list1), L2 = length(list2), common = c1, enrichment = enrichment, pvalue = pvalue, CI = CI))
}


sampleTypes <- c("IDHwt", "IDHmut")
covariates <- "technical"

test = "interactions"; use.bootstraped = T
method = "JADE"; dims = "10"
#test = "components"
#test = "genes"
remove.bcells <- T

fisher.result <- list()
for (sampleType in sampleTypes) {
	
	file.discovery <- glue("ICA-{method}-{dims}.allInteractions-{sampleType}-TCGA-technical.txt")
	file.validation <- glue("ICA-{method}-{dims}.allInteractions-{sampleType}-CGGA693-technical.txt")
	
	tcga <- read.table(file.discovery, sep = "\t", header = T)
	cgga <- read.table(file.validation, sep = "\t", header = T)
	
	tcga <- na.omit(tcga)
	cgga <- na.omit(cgga)
	
	if (remove.bcells == T) {
		indexes <- grep("Bcell", tcga$cellTypes, invert = T)
		tcga <- tcga[indexes, ]
		
		indexes <- grep("Bcell", cgga$cellTypes, invert = T)
		cgga <- cgga[indexes, ]
	}
	universe <- nrow(tcga) * 2 ##approximate
	
	tcga <- tcga[tcga$pvalue < 0.05, ] %>% na.omit()
	cgga <- cgga[cgga$pvalue < 0.05, ] %>% na.omit()
	
	tcga$prognosis <- ifelse(tcga$HR > 0, "Worst", "Better")
	cgga$prognosis <- ifelse(cgga$HR > 0, "Worst", "Better")
	
	tcga.string <- paste(tcga$cellTypes, tcga$components, tcga$prognosis, tcga$bin, sep = "_")
	cgga.string <- paste(cgga$cellTypes, cgga$components, cgga$prognosis, cgga$bin, sep = "_")
	
	fisher.out <- fisherTestFunction(list1 = tcga.string, list2 = cgga.string, universe = universe)
	
	fisher.result[[sampleType]] <- fisher.out
	
	
	if (use.bootstraped == T) {
		file.discovery <- glue("ICA-{method}-{dims}.randomInteractions-{sampleType}-TCGA-technical.txt")
	
		random.data <- read.table(file.discovery, sep = "\t", header = T)
		random.data <- na.omit(random.data)
		
		if (remove.bcells == T) {
			indexes <- grep("Bcell", random.data$cellTypes, invert = T)
			random.data <- random.data[indexes, ]
		}
		
		random.list <- split(random.data, random.data$iteration)
	
		fisher.iteration <- list()
		for (iteration in names(random.list)) {
			random <- random.list[[iteration]]
			universe <- nrow(random) * 2 ##approximate
			
			random <- random[random$pvalue < 0.05, ] %>% na.omit()
			random$prognosis <- ifelse(random$HR > 0, "Worst", "Better")
			random.string <- paste(random$cellTypes, random$components, random$prognosis, random$bin, sep = "_")
			
			fisher.out <- fisherTestFunction(list1 = random.string, list2 = cgga.string, universe = universe)
			fisher.iteration[[iteration]] <- fisher.out
		}				
	}
}

ggTheme <- theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5), #face = "bold"),
	axis.text.y = element_text(size = 12, face = "bold"),
	strip.background =element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 12),
	legend.title = element_blank(), legend.background = element_blank(), legend.position = "none",
	axis.title = element_text(size = 12), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	
	
fisherDf <- bind_rows(fisher.result, .id = "IDH") %>% set_colnames(c("IDH", "TCGA", "CGGA", "Common", "Odds.Ratio", "p.value", "lowerCI", "upperCI"))
fisherDf$Odds.Ratio.log <- log2(fisherDf$Odds.Ratio)
fisherDf$lowerCI.log <- log2(fisherDf$lowerCI)
fisherDf$upperCI.log <- log2(fisherDf$upperCI)

barplot <- ggplot(data = fisherDf, aes(x = IDH, y = Odds.Ratio.log)) + theme_bw() + ylab("Odds Ratio (log2)") + ggTheme +
				geom_bar(stat="identity", fill="white", col = "black") +  scale_x_discrete(label = function(x) stringr::str_wrap(x, width = 10)) + 
				geom_text(aes(label = formatC(p.value, format = "e", digits = 2),  size = 10, y = ifelse(Odds.Ratio.log > 0, Odds.Ratio.log + 0.5, Odds.Ratio.log - 0.5)), color = "black", position = position_dodge(0.9), size = 3.5) +
				geom_errorbar(aes(x = IDH, ymin = lowerCI.log, ymax = upperCI.log), width = 0.1, alpha = 0.9, size =  0.5)
					
svglite::svglite(glue("plots/TCGA.CGGA.oddsRatio.for{test}.svg"), height = 4, width = 2.5)
barplot
dev.off()

				
fisherDf.random <- bind_rows(fisher.iteration, .id = "IDH") %>% set_colnames(c("IDH", "TCGA", "CGGA", "Common", "Odds.Ratio", "p.value", "lowerCI", "upperCI"))
fisherDf.random$Odds.Ratio.log <- log2(fisherDf.random$Odds.Ratio)
fisherDf.random$lowerCI.log <- log2(fisherDf.random$lowerCI)
fisherDf.random$upperCI.log <- log2(fisherDf.random$upperCI)


ggTheme <- theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5), #face = "bold"),
	axis.text.y = element_text(size = 12, face = "bold"),
	strip.background =element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 12),
	legend.title = element_blank(), legend.background = element_blank(), legend.position = "none",
	axis.title = element_text(size = 12), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	

box.plot <- ggplot(data = fisherDf.random, aes(x = 0, y = Odds.Ratio.log)) + geom_boxplot() + geom_point() + theme_bw() + ggTheme +
	coord_cartesian(ylim = c(-2.5, 1.5)) + geom_point(x = 0, y = fisherDf$Odds.Ratio.log, size = 5, col = "red")
	
 
svglite::svglite(glue("plots/TCGA.CGGA.oddsRatio.randomInteractions.svg"), height = 4, width = 2.5)
box.plot
dev.off()


fisherDf.random$IDH <- paste("random", fisherDf.random$IDH, sep = ".")

data.plot <- rbind(fisherDf, fisherDf.random)


ggTheme <- theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1), #face = "bold"),
	axis.text.y = element_text(size = 12, face = "bold"),
	strip.background =element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 12),
	legend.title = element_blank(), legend.background = element_blank(), legend.position = "none",
	axis.title = element_text(size = 12), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	

barplot <- ggplot(data = data.plot, aes(x = IDH, y = Odds.Ratio.log)) + theme_bw() + ylab("Odds Ratio (log2)") + ggTheme +
				geom_bar(stat="identity", fill="white", col = "black") +  scale_x_discrete(label = function(x) stringr::str_wrap(x, width = 10)) + 
				geom_text(aes(label = formatC(p.value, format = "e", digits = 2),  size = 9, y = ifelse(Odds.Ratio.log > 0, Odds.Ratio.log + 0.5, Odds.Ratio.log - 0.5)), color = "darkred", position = position_dodge(0.9), size = 3.5) +
				geom_errorbar(aes(x = IDH, ymin = lowerCI.log, ymax = upperCI.log), width = 0.1, alpha = 0.9, size =  0.5)
					

svglite::svglite(glue("plots/TCGA.CGGA.oddsRatio.actualVsRandom.svg"), height = 4, width = 6)
barplot + ggTheme
dev.off()


############################################
############################################
############################################

############################################Filter Interactions############################################



sampleType <- "mut"; cohort = "TCGA"
dims = 10; method = "JADE"; covariateType = "technical"

file.discovery <- glue("ICA-{method}-{dims}.FDR0.25.robustInteractions-IDH{sampleType}-{cohort}-{covariateType}.txt")
file.validation <- glue("ICA-{method}-{dims}.allInteractions-IDH{sampleType}-CGGA693-technical.txt")


tcga <- read.table(file.discovery, sep = "\t", header = T)
cgga <- read.table(file.validation, sep = "\t", header = T)



tcga <- na.omit(tcga)
cgga <- na.omit(cgga)

tcga$prognosis <- ifelse(tcga$HR > 0, "Worst", "Better")
cgga$prognosis <- ifelse(cgga$HR > 0, "Worst", "Better")

tcga$string <- paste(tcga$cellTypes, tcga$components, tcga$bin, sep = "_")
cgga$string <- paste(cgga$cellTypes, cgga$components, cgga$bin, sep = "_")


tcga.cgga <- merge(tcga, cgga, by = "string")
write.table(tcga.cgga, glue("ICA-{method}-{dims}.FDR0.25.robustInteractions-IDH{sampleType}-tcga.cgga-{covariateType}.txt"), sep = "\t", quote = F, row.names = F)

######1. Within cohort reproducibility filter########
data.plot <- tcga.cgga[tcga.cgga$FDR.x < 0.20 & tcga.cgga$Reproduciblity >= 8, ]

######1. significance in the second cohort filter########
data.plot <- tcga.cgga[tcga.cgga$pvalue.y < 0.1, ]



boxplot(HR.y ~ prognosis.x, data = data.plot, outline = F)

svglite::svglite(glue("plots/TCGA.CGGA.boxplot.for{test}.svg"), height = 4, width = 2.5)
boxplot(HR.y ~ prognosis.x, data = tcga.cgga, outline = F)
dev.off()

###################Other plots######################
data.plot <- tcga.cgga[tcga.cgga$FDR.x < 0.2 & tcga.cgga$Reproduciblity >= 8, ]

svglite::svglite(glue("plots/TCGA.CGGA.boxplot.for{test}-finalInteractions.svg"), height = 4, width = 2.5)
ggplot(data = tcga.cgga, aes(x = prognosis.x, y = HR.y)) + geom_boxplot() + theme_bw() + xlab("Prognosis") + ylab("HR (ln)") + 
	ggpubr::stat_compare_means(label = "p.format", label.x = 1.5)
dev.off()
