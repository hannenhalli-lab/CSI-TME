#library(clusterProfiler)
library(dplyr); library(magrittr);
library(reshape2); library(ggplot2)
library(glue); library(doParallel)
library(igraph); library(svglite)

convert.interactionFormat <- function(interaction.data, input.type = "cellTypePairs") {
	if (input.type == "cellTypePairs") {
		interaction.data <- interaction.data %>% tidyr::separate(., col = cellTypes, sep = ":", into = c("cellType1", "cellType2")) %>% 
				tidyr::separate(., col = components, sep = ":", into = c("component1", "component2"))
	}
	return(interaction.data)
}


idh = "mut"
cohort <- "TCGA"
dims = "10"


interaction.file <- glue("ICA-JADE-{dims}.FDR0.25.robustInteractions-IDH{idh}-tcga.cgga-technical.txt")
interaction.data <- read.table(interaction.file, sep = "\t", header = T)
names(interaction.data)
interaction.data$prognosis <- ifelse(interaction.data$HR.x > 0, "Worse", "Better")
interaction.data <- interaction.data %>% subset(., FDR.x < 0.20 & Reproduciblity > 7)


interaction.data %>% split(.,.$prognosis.x) %>% lapply(., function(x) quantile(x$HR.x))

Interaction.Count <- interaction.data %>% group_by(prognosis) %>% tally()

svglite("plots/barplot-numberOfinteractions-TCGA.svg", height = 2.5, width = 1.5)
ggplot(data = Interaction.Count, aes(x = prognosis, y = n)) + 
		geom_bar(position="dodge", stat="identity", fill = "white", col = "black") + theme_bw()
dev.off()		



bin.distribution <- interaction.data %>% group_by(prognosis, bin.x) %>% tally()
total.n <- nrow(interaction.data)
total.prognosis.n <- table(interaction.data$prognosis) %>% data.frame()
bin.distribution$proportion <- bin.distribution$n / total.n


bin.plot.data <- reshape2::dcast(bin.distribution, prognosis ~ bin.x, value.var = 'n')
rownames(bin.plot.data) <- bin.plot.data[ ,1]
bin.plot.data <- bin.plot.data[ ,-1]
rownames(bin.plot.data) <- c("Anti-tumor", "Pro-tumor")

svglite::svglite(glue("plots/ICA-JADE-{dims}.FDR0.20.robustInteractionsinteractions.bin.distribution.svg"), height = 1, width = 2.5) 
pheatmap::pheatmap(bin.plot.data, display_numbers = bin.plot.data, cluster_rows = F, cluster_cols = F, fontsize_number = 16)
dev.off()


bin.distribution <- merge(bin.distribution, total.prognosis.n, by.x = "prognosis", by.y = "Var1")
bin.distribution$fraction <- bin.distribution$n / bin.distribution$Freq
bin.plot.data <- reshape2::dcast(bin.distribution, prognosis ~ bin.x, value.var = 'fraction')
rownames(bin.plot.data) <- bin.plot.data[ ,1]
bin.plot.data <- bin.plot.data[ ,-1]
rownames(bin.plot.data) <- c("Anti-tumor", "Pro-tumor")




bin.plot.data <- round(bin.plot.data, 2)
svglite::svglite(glue("plots/ICA-JADE-{dims}.FDR0.20.robustInteractionsinteractions.bin.distribution-fraction.svg"), height = 1, width = 2.5) 
pheatmap::pheatmap(bin.plot.data, display_numbers = bin.plot.data, cluster_rows = F, cluster_cols = F, fontsize_number = 16)
dev.off()




