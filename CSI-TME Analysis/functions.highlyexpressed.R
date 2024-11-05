library(magrittr); library(glue); library(dplyr)
library(clusterProfiler); library(org.Hs.eg.db)
library(preprocessCore)
orgDb <- "org.Hs.eg.db"

source('../../melanomaSplicing/generalFunctions.R')


expression.path <- "tcga_mut_filtered/"
files <- list.files(expression.path)

gene.expression.list <- list()
for (file.name in files) {
	gene.expression <- read.table(glue("{expression.path}/{file.name}"), header = T, row.names = 1)
	gene.expression <- log2(gene.expression + 0.001)
	gene.expression <- normalize.quantiles(data.matrix(gene.expression)) %>% 
		set_rownames(rownames(gene.expression)) %>% set_colnames(colnames(gene.expression))
	
	cellType <- gsub("\\.txt", "", file.name)
	gene.expression.list[[cellType]] <- gene.expression
	message(glue("read data for {file.name}"))
}

gene.median.list <- list()
for (cellType in names(gene.expression.list)) {
	gene.median <- apply(gene.expression.list[[cellType]], 1, median, na.rm = T)
	gene.median <- gene.median[order(gene.median, decreasing = T)]
	gene.median.list[[cellType]] <- gene.median
}


#top.genes.n <- 1000
#gene.list <- lapply(gene.median.list, function(x) x[1:top.genes.n] %>% names())

top.genes.p <- 10
gene.list <- lapply(gene.median.list, function(x) {size.x = length(x); top.genes.n = floor((size.x)* top.genes.p) / 100; x[1:top.genes.n] %>% names()})


gene.frequency <- unlist(gene.list) %>% table
gene.specific <- gene.frequency[gene.frequency < 3] %>% names()
gene.ubiqutous <- gene.frequency[gene.frequency >= 3] %>% names()
gene.list <- lapply(gene.list, function(x) x[x %in% gene.specific])


#background.genes <- lapply(gene.median.list, names) %>% unlist %>% unique()

#common.genes <- Reduce(intersect, gene.list)

enrichment.list <- lapply(gene.list, function(x) {
	if (length(x) < 3) {
		functionalEnrichment("GAPDH", gene.id = "SYMBOL", background = NULL, T2G = NULL, T2N = NULL, method = "GO", ontology = "BP", method.adjust = 'fdr', organism = "human", p = 0.05, q = 1)
	} else {
		functionalEnrichment(x, gene.id = "SYMBOL", background = NULL, T2G = NULL, T2N = NULL, method = "GO", ontology = "BP", method.adjust = 'fdr', organism = "human", p = 0.05, q = 1)
	}
})
enrichment.list.simplified <- lapply(enrichment.list, function(x) simplify(x, cutoff = 0.7))


enrichment <- enrichment.list.simplified %>% lapply(data.frame) %>% bind_rows(.id = "cellType")
enrichment <- data.frame(enrichment) %>% mutate(., cellType = cellType, .before = 1)

enrichment <- enrichment %>% subset(., geneID != "GAPDH" & p.adjust < 0.20)
#enrichment$signif <- -log10(enrichment$qvalue)
enrichment$signif <- -log10(enrichment$p.adjust)

foreground <- lapply(enrichment$GeneRatio, function(x) eval(parse(text = x))) %>% unlist()
background <- lapply(enrichment$BgRatio, function(x) eval(parse(text = x))) %>% unlist()

enrichment$OE <- log2(foreground / background)

enrichment.subset <- enrichment %>%  arrange(desc(signif)) %>% group_by(cellType) %>% dplyr::slice(1:20)

#enrichment.subset <- enrichment %>%  arrange(desc(OE)) %>% group_by(Dimension) %>% dplyr::slice(1:7)
					  
data.plot <- reshape2::dcast(data = enrichment.subset,formula = Description ~ cellType, value.var = "OE")
data.plot <- data.plot[ ,-1] %>% set_rownames(data.plot[ ,1])
data.plot[is.na(data.plot)] <- 0
col.sum <- colSums(data.plot)
data.plot <- data.plot[ ,col.sum > 0]

pheatmap::pheatmap(data.plot, scale = 'none', cell.width = 2, cell.height = 2, clustering_method = "ward.D2", fontsize_row = 6, show_rownames = T)


output.file <- glue("plots/heatmap.GO.enrichment-highlyExpressed-multiTest.cellTypes-TCGA-IDHmut.svg")
svglite::svglite(file = output.file, height = 8, width = 6)
pheatmap::pheatmap(data.plot, scale = 'none', cell.width = 2, cell.height = 2, clustering_method = "ward.D2", fontsize_row = 6, show_rownames = T)
dev.off()



######################################################################################################
gene.list <- lapply(gene.median.list, function(x) {size.x = length(x); top.genes.n = floor((size.x)* top.genes.p) / 100; x[1:top.genes.n] %>% names()})
gene.list <- lapply(gene.list, function(x) data.frame(gene = x))
gene.data <- bind_rows(gene.list, .id = "cellType.expression")


surv.file <- "survivalAnalysis.Genes-IDHmut-TCGA-cellTypes.txt"
surv.data <- read.table(surv.file, sep = "\t", header = T)

data.plot <- merge(surv.data, gene.data, by.x = "Component", by.y = "gene")

if (prognosis == "positive") {
	surv.genes <- surv.data[surv.data$HR != 0 & surv.data$FDR < 0.05, ]
	gene.list <- split(surv.genes, surv.genes$cellType) %>% lapply(., function(x) x$Component)
	
}
if (prognosis == "negative") {
	surv.genes <- surv.data[surv.data$HR < 0 & surv.data$FDR < 0.05, ]
	gene.list <- split(surv.genes, surv.genes$cellType) %>% lapply(., function(x) x$Component)
	
}


