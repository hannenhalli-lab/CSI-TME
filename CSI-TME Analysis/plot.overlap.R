###########plot#############
library(dplyr); library(magrittr); library(glue)
library(doParallel); library(tidyr); library(ggplot2)
source('../functionsForInteractionAnalysis.R')

ggTheme <- theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1), axis.text.y = element_text(size = 10),
			axis.title = element_text(size = 10, face ="bold"), panel.border = element_rect(color = "grey", fill = NA, size = 2),
			panel.grid.major = element_blank(), #panel.grid.minor.x = element_line(colour="black", size=1)
			strip.text.x = element_text(size = 10, colour = "brown")
			)



tooMany <- F
input.file <- "enrichmentValues.IFN.ICs.negativeGenes.Tcells.txt"
output.file <- "plots/dotplot.overlap.IFN.negativeGenes.Tcell.svg"


overlap.df <- read.table(input.file, sep = "\t", header = T)



overlap.df$significance <- ifelse(overlap.df$FDR < 0.20 & overlap.df$enrichment > 1, "Significant", "Not-Significant")
overlap.df$component <- factor(overlap.df$component);
states <- overlap.df[overlap.df$FDR < 0.20, 'State'] %>% unlist() %>% unique()

overlap.df$IC <- gsub("Dim", "IC", overlap.df$component)
overlap.df$IC <- factor(overlap.df$IC, levels = unique(overlap.df$IC))

if (tooMany) {
	overlap.df <- overlap.df[overlap.df$State %in% states, ]
}


overlap.df$xLines <- as.numeric(factor(overlap.df$IC)) + 0.5
overlap.df$yLines <- as.numeric(factor(overlap.df$State)) + 0.5
overlap.df$log.enrich <- log2(overlap.df$enrichment + 0.01)


#breaks = c(-6, -4, -2, 0, 2, 4, 6)
#breaks = c(0, 2, 5, 10, 100)
box.plot <- ggplot(overlap.df, aes(x = IC, y = State)) + theme_bw() + ggTheme +
	 		geom_point(aes(size = enrichment, shape = significance), colour = "black") + 
			scale_shape_manual(values=c(21, 16)) + scale_size_continuous(range = c(0,10), breaks = c(0, 2, 5, 10, 100)) +
			geom_hline(aes(yintercept = yLines), color = "grey") + geom_vline(aes(xintercept = xLines), color = "grey")
 
svglite::svglite(file = output.file, height = 2, width = 6)
box.plot
dev.off()