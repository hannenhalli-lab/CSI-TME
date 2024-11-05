library(dplyr); library(magrittr); library(glue)
library(doParallel); library(tidyr)
source('../functionsForInteractionAnalysis.R')

geneset.file <- "geneSets.Compilation.xlsx"
idh = "mut"; target.cellType <- "Tcell"
cohort <- "TCGA"
dims = "10"
load(glue("ICA-JADE-{dims}.signatureGenes.{tolower(cohort)}.IDH{idh}--filtered-noNorm-optimalRank.Rda"))
load(glue("ICA-JADE-{dims}.cellTypes.{tolower(cohort)}.IDH{idh}-filtered-noNorm-optimalRank.Rda"))


universe.list <- list()
for(cellType in names(ica.res.list)) {
	universe <- ica.res.list[[cellType]]$S %>% rownames()
	universe.list[[cellType]] <- universe
}
signatures.cellType <- signature.list[[target.cellType]]
universe <- universe.list[[target.cellType]]

sheets <- readxl::excel_sheets(geneset.file)
#sheets <- c("Venteicher")
#sheets <- c("Wechter.TimeCourse")
#sheets <- c("Tcell.Proliferation");
sheets <- "Dual.Tcell"

#sheets <- c("Biomarkers");

#sheets <-c("SnG.Clusters");
#sheets <-c("IDH1.correlation");
#sheets <- c("Telomeres")
#sheets <- "lm22"
#sheets <- c("Wechter")
#sheets <- c("CellAge")
gene.sets <- lapply(sheets, function(x) readxl::read_excel(geneset.file, sheet = x))
names(gene.sets) <- sheets


fisher.result.list <- list(); ic.genes <- "negative"
for (geneset in names(gene.sets)) {
	gene.list <- gene.sets[[geneset]] %>% split(., .$State) %>% lapply(., function(x) x$Gene)
	
	res.geneset.list <- list()
	for (set.name in names(gene.list)) {
		set.genes <- gene.list[[set.name]] %>% toupper
		list.res.comps <- list()
		for (component in names(signatures.cellType)) {
			component.genes <- signatures.cellType[[component]][[ic.genes]]
			if (length(component.genes) > 0) {
				fisher.result <- fisherTestFunction(component.genes, set.genes, universe = universe)
				list.res.comps[[component]] <- fisher.result
			}
		}
		res <- bind_rows(list.res.comps, .id = "component")
		res.geneset.list[[set.name]] <- res
	}
	fisher.df <- bind_rows(res.geneset.list, .id = "State")
	fisher.df$FDR <- p.adjust(fisher.df$pvalue, method = "fdr")
	fisher.result.list[[geneset]] <- fisher.df
}

fisher.result.df <- bind_rows(fisher.result.list, .id = "geneset")
fisher.result.df[fisher.result.df$FDR < 0.2, ]
#fisher.result.df[fisher.result.df$State == "CD8_c1_Tex", ] %>% data.frame
write.table(fisher.result.df, file = glue("enrichmentValues.{sheets}.ICs.{target.cellType}.txt"), sep = "\t", quote = F, row.names = F)


file.name <- glue("fisher.enrichment.genesets.ICA-JADE-{dims}.signatureGenes.{tolower(cohort)}.IDH{idh}-{target.cellType}-{sheets}.xlsx")
if (!file.exists(file.name)) {
	xlsx::write.xlsx(fisher.result.df, file.name, sheetName = ic.genes)
} else {
	xlsx::write.xlsx(fisher.result.df, file.name, sheetName = ic.genes, append = TRUE) 
}


###########plot#############
library(dplyr); library(magrittr); library(glue)
library(doParallel); library(tidyr)
source('../functionsForInteractionAnalysis.R')

library(ggplot2)
idh = "mut";
cohort <- "TCGA"
dims = "10"
target.cellType <- "Tcell"; 

file.name <- glue("fisher.enrichment.genesets.ICA-JADE-{dims}.signatureGenes.{tolower(cohort)}.IDH{idh}-{target.cellType}.xlsx")
sheets <- readxl::excel_sheets(file.name)
overlap.list <- lapply(sheets, function(x) readxl::read_excel(file.name, sheet = x))
names(overlap.list) <- sheets

ggTheme <- theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1), axis.text.y = element_text(size = 10),
			axis.title = element_text(size = 10, face ="bold"), panel.border = element_rect(color = "grey", fill = NA, size = 2),
			panel.grid.major = element_blank(), #panel.grid.minor.x = element_line(colour="black", size=1)
			strip.text.x = element_text(size = 10, colour = "brown")
			)
######positive######

target.genes <- "negative"; tooMany = F
output.file <- glue("plots/dotplot.enrichment.genesets-{target.genes}.{sheets}.ICA-JADE-{dims}.signatureGenes.{tolower(cohort)}.IDH{idh}-{target.cellType}.svg")

overlap.df <- overlap.list[[target.genes]]
overlap.df$IC <- gsub("Dim", "IC", overlap.df$component)
overlap.df$IC <- factor(overlap.df$IC, levels = unique(overlap.df$IC))
colnames(overlap.df)[7] <- "enrichment"
overlap.df$significance <- ifelse(overlap.df$FDR < 0.20 & overlap.df$enrichment > 1, "Significant", "Not-Significant")
overlap.df$component <- factor(overlap.df$component); overlap.df$geneset <- factor(overlap.df$geneset)
states <- overlap.df[overlap.df$FDR < 0.20, 'State'] %>% unlist() %>% unique()

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
 
svglite::svglite(file = output.file, height = 2.5, width = 8)
box.plot
dev.off()


########################################################
########################################################
########################negative########################


target.genes <- "negative"
output.file <- glue("plots/dotplot.enrichment.genesets-{target.genes}.ICA-JADE-{dims}.signatureGenes.{tolower(cohort)}.IDH{idh}-{target.cellType}.svg")


overlap.df <- overlap.list[['negative']] %>% subset(., geneset != "Macrophages")
overlap.df$significance <- ifelse(overlap.df$FDR < 0.25 & overlap.df$enrichment > 1, "Significant", "Not-Significant")
overlap.df$component <- factor(overlap.df$component); overlap.df$geneset <- factor(overlap.df$geneset)
overlap.df$xLines <- as.numeric(factor(overlap.df$component)) + 0.5
overlap.df$yLines <- as.numeric(factor(overlap.df$State)) + 0.5

box.plot <- ggplot(overlap.df, aes(x = component, y = State, color = enrichment)) + theme_bw() + ggTheme +
	 		geom_point(aes(size = enrichment, shape = significance), colour = "black") + 
			scale_shape_manual(values=c(21, 16)) + scale_size_continuous(range = c(1,10), breaks = c(0, 2, 5, 10, 100)) +
			geom_hline(aes(yintercept = yLines), color = "grey") + geom_vline(aes(xintercept = xLines), color = "grey")
 


svg(file = output.file, height = 9, width = 6)
box.plot
dev.off()