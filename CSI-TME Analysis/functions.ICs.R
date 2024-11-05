library(clusterProfiler);
library(dplyr); library(magrittr);
library(reshape2); library(ggplot2);
library(glue); library(NMF)
library(org.Hs.eg.db)
orgDb <- "org.Hs.eg.db"

source('../../melanomaSplicing/generalFunctions.R')



sampleType = "mut"
cohort <- "TCGA"; dims = "10"
load(glue("ICA-JADE-{dims}.signatureGenes.tcga.IDHmut--filtered-noNorm-optimalRank.Rda"))


#gene.indexes <-  lapply(nmf.optimal.list, extractFeatures)
#gene.sets <- lapply(nmf.optimal.list, genesets.nmf)

set <- "msigdb"; category = "H"; genes = "negative"
#set <- "GO"; genes = "positive"
geneSets = msigdbr::msigdbr(species = "human", category = category)
term2Gene = geneSets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

enrichment.list <- list()
for(cellType in names(signature.list)) {
	if (genes == "positive") {
		gene.list <- signature.list[[cellType]] %>% lapply(., function(x) x$positive)
		
	}
	if (genes == "negative") {
		gene.list <- signature.list[[cellType]] %>% lapply(., function(x) x$negative)
	}
	
	
	
	if (set == "kegg") {
		enrichment <- lapply(gene.list, function(x) {
			if (length(x) < 3) {
				x <- "GAPDH"
				x <- mapIds(org.Hs.eg.db, x, 'ENTREZID', 'SYMBOL')
				functionalEnrichment(x, gene.id = "ncbi-geneid", background = NULL, T2G = NULL, T2N = NULL, method = "KEGG", method.adjust = 'none', organism = "human", p = 0.05, q = 0.2)
			} else {
				x <- mapIds(org.Hs.eg.db, x, 'ENTREZID', 'SYMBOL')
				functionalEnrichment(x, gene.id = "ncbi-geneid", background = NULL, T2G = NULL, T2N = NULL, method = "KEGG", method.adjust = 'none', organism = "human", p = 0.05, q = 0.2)
			}
		})
	}
	if (set == "GO") {
		enrichment <- lapply(gene.list, function(x) {
			if (length(x) < 3) {
				functionalEnrichment("GAPDH", gene.id = "SYMBOL", background = NULL, T2G = NULL, T2N = NULL, method = "GO", ontology = "BP", method.adjust = 'fdr', organism = "human", p = 0.05, q = 1)
			} else {
				functionalEnrichment(x, gene.id = "SYMBOL", background = NULL, T2G = NULL, T2N = NULL, method = "GO", ontology = "BP", method.adjust = 'fdr', organism = "human", p = 0.05, q = 1)
			}
		})
	}
	if(set == "msigdb") {
		enrichment <- lapply(gene.list, function(x) {
			if (length(x) < 3) {
				functionalEnrichment("GAPDH", gene.id = "SYMBOL", background = NULL, T2G = term2Gene, T2N = NULL, method = "General",  method.adjust = 'fdr', organism = "human", p = 0.05, q = 0.2)
			} else {
				functionalEnrichment(x, gene.id = "SYMBOL", background = NULL, T2G = term2Gene, T2N = NULL, method = "General",  method.adjust = 'fdr', organism = "human", p = 0.05, q = 0.2)
			}
		})
	}
	
	enrichment.list[[cellType]] <- enrichment
	print(glue("done for {cellType}"))
}

if (set == "GO") {
	save(enrichment.list, file = glue("{set}.enrichment.ICs-{genes}.{dims}-{cohort}-IDH{sampleType}.Rda"))
}
if (set == "msigdb") {
	save(enrichment.list, file = glue("{set}-{category}.enrichment.ICs-{genes}.{dims}-{cohort}-IDH{sampleType}.Rda"))
}

library(pheatmap)
#####load(glue("GO.enrichment.ICs.{dims}-{cohort}-IDH{sampleType}.Rda"))
plot.list <- list()
for(cellType in names(enrichment.list)) {
	if (set == "GO") {
		enrichment <- enrichment.list[[cellType]] %>% lapply(., function(x) simplify(x, cutoff = 0.7)) %>% lapply(data.frame) %>% bind_rows(.id = "Dimension")
	} else {
		enrichment <- enrichment.list[[cellType]] %>% lapply(data.frame) %>% bind_rows(.id = "Dimension")
	}
	
	enrichment <- data.frame(enrichment) %>% mutate(., cellType = cellType, .before = 1)
	enrichment <- enrichment %>% subset(., geneID != "GAPDH" & p.adjust < 0.2)
	#enrichment$signif <- -log10(enrichment$qvalue)
	enrichment$signif <- -log10(enrichment$p.adjust)
	
	foreground <- lapply(enrichment$GeneRatio, function(x) eval(parse(text = x))) %>% unlist()
	background <- lapply(enrichment$BgRatio, function(x) eval(parse(text = x))) %>% unlist()
	
	enrichment$OE <- log2(foreground / background)
	
	enrichment.subset <- enrichment %>%  arrange(desc(signif)) %>% group_by(Dimension) %>% dplyr::slice(1:20)

	#enrichment.subset <- enrichment %>%  arrange(desc(OE)) %>% group_by(Dimension) %>% dplyr::slice(1:7)
						  
  	data.plot <- dcast(data = enrichment.subset,formula = Description ~ Dimension, value.var = "OE")
  	data.plot <- data.plot[ ,-1] %>% set_rownames(data.plot[ ,1])
	data.plot[is.na(data.plot)] <- 0
	col.sum <- colSums(data.plot)
	data.plot <- data.plot[ ,col.sum > 0]
  	plot.list[[cellType]] <- pheatmap(data.plot, scale = 'none', cell.width = 2, cell.height = 2, clustering_method = "ward.D2", fontsize_row = 6, show_rownames = T)
	
	message(glue("done for {cellType}"))
}

save(plot.list, file = glue("listOfFunctionalHeatmaps-{set}.enrichment-multiTest.ICs-{genes}.{dims}-{cohort}-IDH{sampleType}.Rda"))




save_pheatmap_svg <- function(x, filename, width = 5, height = 6) {
	####modified from https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   svglite::svglite(filename, width = width, height = height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}


cellType <- "Tcell"
output.file <- glue("plots/heatmap.{set}.enrichment-multiTest.ICs-{cellType}.{genes}.{dims}-{cohort}-IDH{sampleType}.svg")

save_pheatmap_svg(plot.list[['Tcell']], output.file)

########################################################################







plot.list <- list()
for(cellType in names(enrichment.list)) {
	enrichment <- enrichment.list[[cellType]] %>% lapply(data.frame) %>% bind_rows(.id = "Dimension")
	
	enrichment <- data.frame(enrichment) %>% mutate(., cellType = cellType, .before = 1)
	enrichment <- enrichment %>% subset(., geneID != "GAPDH" & p.adjust < 0.2)
	enrichment$signif <- -log10(enrichment$p.adjust)
	
	foreground <- lapply(enrichment$GeneRatio, function(x) eval(parse(text = x))) %>% unlist()
	background <- lapply(enrichment$BgRatio, function(x) eval(parse(text = x))) %>% unlist()
	
	enrichment$OE <- log2(foreground / background)
	
	
	
	data.plot <- dcast(data = enrichment,formula = Description ~ Dimension, value.var = "OE")
	data.plot <- data.plot[ ,-1] %>% set_rownames(data.plot[ ,1])
	data.plot[is.na(data.plot)] <- 0
	
	plot.list[[cellType]] <- pheatmap(data.plot, scale = 'none', cell.width = 2, cell.height = 2, clustering_method = "ward.D2", fontsize_row = 6, show_rownames = T)
	#enrichmentDf.list[[cellType]] <- enrichment
	#enrichment.Df <- rbind(enrichment.Df, enrichment)
}

save(plot.list, file = glue("listOfFunctionalHeatmaps-{set}-{category}.enrichment-multiTest.ICs-{genes}.{dims}-{cohort}-IDH{sampleType}.Rda"))

output.file <- glue("plots/heatmap.{set}.enrichment-multiTest.ICs-{genes}.{dims}-{cohort}-IDH{sampleType}.svg")
svglite::svglite(file = output.file, height = 6, width = 6)
plot.list[['Tcell']]
dev.off()  


