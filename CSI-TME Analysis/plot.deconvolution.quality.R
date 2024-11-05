library(magrittr); library(glue); library(dplyr)

source('../../melanomaSplicing/generalFunctions.R')
source('../functionsForInteractionAnalysis.R')

overlap.file <- "overlap.deconvovledTopGenes.cellTypes.txt"
survival.gene.file <- "correlation.survival.genes.cellTypes.txt"
survival.component.file <- "correlation.survival.IC.cellTypes.txt"


overlap <- read.table(overlap.file, sep = "\t", header = T)
survival.gene <- read.table(survival.gene.file, sep = "\t", header = T)
survival.component <- read.table(survival.component.file, sep = "\t", header = T)

overlap.data.plot <- overlap[overlap$type == "same cell type", ]
survival.gene.data.plot <- survival.gene[survival.gene$type == "same cell type", ]
survival.component.data.plot <- survival.component[survival.component$comparison == "same cell type", ]


output.file <- glue("plots/dotplot.overlap.topGenes.TCGA.CGGA-SameCellType.svg")
svglite::svglite(file = output.file, height = 3.5, width = 2.5)
ggplot(data = overlap.data.plot, aes(x = "overlap", y = TCGA)) + geom_point(aes(size=oddsRatio)) + scale_size_continuous(range = c(2,5), breaks = c(2,8,14)) + theme_bw()
dev.off()



#data.heatmap <- merge(survival.gene.data.plot, survival.component.data.plot, by.x = "cellType1", by.y = "cellType")
#data.heatmap <- data.heatmap[ ,c("correlation.x", "correlation.y")] %>% set_rownames(data.heatmap$cellType1) %>% set_colnames(c("Genes", "ICs"))
#pheatmap::pheatmap(data.heatmap, scale = 'none', cluster_rows = F, cluster_cols = F)

data.plot1 <- survival.gene.data.plot[ ,c("cellType1", "correlation")] %>% mutate(feature = "Gene") %>% set_names(c("cellType", "correlation", "feature"))
data.plot2 <- survival.component.data.plot[ ,c("cellType", "correlation")] %>% mutate(feature = "ICs") %>% set_names(c("cellType", "correlation", "feature"))

data.plot <- rbind(data.plot1, data.plot2)


output.file <- glue("plots/dotplot.correlation.survival.TCGA.CGGA-SameCellType.svg")
svglite::svglite(file = output.file, height = 3.5, width = 3.0)
ggplot(data = data.plot, aes(x = feature, y = cellType)) + geom_point(aes(size=correlation)) + scale_size_continuous(range = c(1,6), breaks = c(0.0, 0.2,0.4,0.6, 0.8, 1.0)) + theme_bw()
dev.off()
