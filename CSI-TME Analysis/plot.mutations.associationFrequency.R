			#####################frequency of genes with association################
			

library(dplyr); library(magrittr); library(glue)
library(tidyr); library(ggplot2); library(tidyverse)
source('../functionsForInteractionAnalysis.R')


ggTheme <- theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1), axis.text.y = element_text(size = 10),
			axis.title = element_text(size = 10, face ="bold"), panel.border = element_rect(color = "grey", fill = NA, size = 2),
			panel.grid.major = element_blank(), #panel.grid.minor.x = element_line(colour="black", size=1)
			strip.text.x = element_text(size = 10, colour = "brown")
			)

idh = "mut"
cohort <- "tcga"
dims = 10; method = "JADE"
input.type = "interactions"


interactionFile <- glue("ICA-JADE-{dims}.FDR0.25.robustInteractions-IDHmut-tcga.cgga-technical.txt")
interactions <- read.table(interactionFile, sep = "\t", header = T)
interactions.robust <- interactions %>% subset(., FDR.x < 0.20 & Reproduciblity > 7)
interactions.robust.list <- split(interactions.robust, interactions.robust$prognosis.x) %>% lapply(., function(x) x$string)

if (input.type == "interactions") {
	file.input <- glue("ICA-JADE-{dims}.FDR0.20.robustInteractions-IDH{idh}-association.Mutations.txt")
}
if (input.type == "cellTypes") {
	file.input <- glue("ICA-{method}-{dims}.cellTypes.{cohort}.IDH{idh}-assocaition.Mutations.txt")
}

input.data <- read.table(file.input, sep = "\t", header = T)
input.data$Association <- ifelse(input.data$Estimate > 0 , "Positive", "Negative")

if (input.type == "cellTypes") {
	association.significant <- input.data[input.data$FDR < 0.20,  ]
	
	#####which genes are frequently assocaited across cell types#####
	output.file <-  glue("plots/ICA-JADE-{dims}.association.IC.Mutations-GeneFrequency-barplot.svg")
	
	gene.freq <- table(association.significant$response) %>% data.frame()
	gg.p <- ggplot(data = gene.freq, aes(x = reorder(Var1, -Freq), y = Freq)) + theme_bw() + ggTheme+ 
			geom_bar(stat = "identity", col = "black", fill = "white") + xlab("Gene") + ylab("# Associations")
	
	
	ggsave(file = output.file, height = 3, width = 4)
	gg.p 
	dev.off()
	
			
	#####which cell types are frequently assocaited across cell types#####
	cellType.predictor <- paste(association.significant$cellType, association.significant$predictor, sep = ":")
	#cellType.freq <- table(association.significant$cellType) %>% data.frame()
	#ggplot(data = gene.freq, aes(x = reorder(Var1, -Freq), y = Freq)) + theme_bw() + ggTheme+ 
	#geom_bar(stat = "identity", col = "black", fill = "white") + xlab("Gene") + ylab("# Associations")
	
	association.significant$cellType.predictor <- cellType.predictor
	
	
	
	data.plot <- reshape2::dcast(association.significant, formula = cellType.predictor + Association ~ response, fun.aggregate = NULL, value.var = "significance")
	#data.plot[is.na(data.plot)] <- 0 ###assign 0 if association was not detected
	
	data.plot <- reshape2::melt(data.plot)
	
	
	ggTheme <- theme(axis.text.x = element_text(size = 10, angle = 35, vjust = 1, hjust = 1), axis.text.y = element_text(size = 10),
				axis.title = element_text(size = 10, face ="bold"), panel.border = element_rect(color = "grey", fill = NA, size = 2),
				#panel.grid.major = element_blank(), panel.grid.minor.x = element_line(colour="black", size=1),
				strip.text.x = element_text(size = 10, colour = "brown")
				)
	
	
	output.file <-  glue("plots/ICA-JADE-{dims}.association.IC.Mutations-dotplot.svg")
	
	
	dotplot <- ggplot(data.plot, aes(x = variable, y = cellType.predictor)) + theme_bw() + ggTheme +
		geom_point(aes(size = value)) + facet_wrap(~Association) + xlab("Mutated genes") + ylab("cell type ICs") +
		scale_size_continuous(range = c(1,4), breaks = c(1, 5, 10, 15)) 
							
	ggsave(file = output.file, height = 8, width = 6)
	dotplot 
	dev.off()
	
	
	#####which cell types are frequently assocaited across cell types#####
	cellType.predictor <- paste(association.significant$cellType, association.significant$predictor, sep = ":")
	cellType.predictor <- table(cellType.predictor) %>% data.frame()
	ggplot(data = cellType.predictor, aes(x = reorder(cellType.predictor, -Freq), y = Freq)) + theme_bw() + ggTheme+ 
	 geom_bar(stat = "identity", col = "black", fill = "white") + xlab("cell type ICs") + ylab("# Associations")
	
	
} else {
	
	association.significant <- input.data[input.data$FDR < 0.20 & input.data$Estimate > 0,  ]
	association.significant$significance <- -log10(association.significant$FDR)
	association.significant.list <- split(association.significant, association.significant$response)
	all.genes <- association.significant$response %>% unique()

	universe.size <- nrow(interactions.robust)
	enrichment.list <- list()
	for (prognosis in names(interactions.robust.list)) {
		list1 <- interactions.robust.list[[prognosis]]
		fisher.result.list <- list()
		for (candidate in all.genes) {
			list2 <- association.significant.list[[candidate]]$predictor
			fisher.result <- fisherTestFunction(list1, list2, universe = universe.size)
			fisher.result.list[[candidate]] <- fisher.result
		}
		enrichment.list[[prognosis]] <- bind_rows(fisher.result.list, .id = "Mutated Gene")
	}

	enrichmentDf <- bind_rows(enrichment.list, .id = "Prognosis") %>% data.frame()
	enrichmentDf$FDR <- p.adjust(enrichmentDf$pvalue, method = "fdr")
	enrichmentDf$significance <- ifelse(enrichmentDf$FDR < 0.25, "Significant", "Insignificant")
	inf.index <- !is.finite(enrichmentDf$enrichment.odds.ratio)
	enrichmentDf$enrichment.odds.ratio[inf.index] <- 15 ###assiign max value 
	enrichmentDf$enrichment <- ifelse(enrichmentDf$enrichment.odds.ratio > 1, "Enriched", "Depleted")

	#enrichmentDf <- enrichmentDf[enrichmentDf$FDR < 0.25, ]
	enrichmentDf$Prognosis <- ifelse(enrichmentDf$Prognosis == "Better", "Anti Tumor", "Pro Tumor")



	dotplot <- ggplot(enrichmentDf, aes(x = Prognosis, y = Mutated.Gene)) + theme_bw() + #ggTheme +
				 		geom_point(aes(size = enrichment.odds.ratio, color = enrichment, shape = significance)) +
						#scale_color_gradient2(midpoint = midpoint.value, low = "purple", mid = "darkgreen", high = "red", space = "Lab") +
						scale_shape_manual(values=c(21, 16)) + scale_size_continuous(range = c(1,4), breaks = c(1, 5, 10, 15))  +
						xlab("Interactions") + ylab("Mutated Genes")




	ggTheme <- theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1), axis.text.y = element_text(size = 10),
				axis.title = element_text(size = 10, face ="bold"), panel.border = element_rect(color = "grey", fill = NA, size = 2),
				panel.grid.major = element_blank(), #panel.grid.minor.x = element_line(colour="black", size=1)
				strip.text.x = element_text(size = 10, colour = "brown")
				)
			
	output.file <-  glue("plots/ICA-JADE-{dims}.FDR0.20.robustInteractions-IDH{idh}-association.Mutations-Prognosis-plot.svg")
	svg(file = output.file, height = 7, width = 3.5)
	dotplot + ggTheme
	dev.off()	
}


#################obsolete#######################


data.gene.frequency <- association.significant %>% group_by(response, prognosis) %>% tally() %>%  reshape2::dcast(response ~ prognosis, value.var = 'n')
data.gene.frequency[is.na(data.gene.frequency)] <- 0
#data.gene.frequency <- reshape2::melt(data.gene.frequency)

#ggplot(data = data.gene.frequency, aes(x = `Anti Tumor`, y = `Pro Tumor`)) + geom_point() + geom_text(aes(label = response))
#ggplot(data = data.gene.frequency, aes(x = variable, y = value)) + geom_boxplot() + geom_text(aes(label = response))


data.gene.frequency <- reshape2::melt(data.gene.frequency)
data.gene.frequency <- merge(data.gene.frequency, background, by.x = "variable", by.y = "prognosis")
data.gene.frequency$frequency <- data.gene.frequency$value / data.gene.frequency$n

data.gene.frequency <- data.gene.frequency %>% reshape2::dcast(response ~ variable, value.var = 'frequency')
ggplot(data = data.gene.frequency, aes(x = `Anti Tumor`, y = `Pro Tumor`)) + geom_point() + geom_text(aes(label = response))
