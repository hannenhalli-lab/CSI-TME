library(glue); library(dplyr);
library(magrittr); library(tximport);
library(svglite); library(GenomicRanges )
library(gUtils); library(ggplot2); 
library(ggpubr); library(gridExtra)
#load('NMF.signatureGenes.tcga.IDHmut--filtered-noNorm-optimalRank.Rda')

idh = "mut"
cohort <- "TCGA"
dims = "10"
load(glue("ICA-JADE-{dims}.signatureGenes.{tolower(cohort)}.IDH{idh}--filtered-noNorm-optimalRank.Rda"))


sample.genes <- function(chromosomes, gff.gr) {
	lapply(names(chromosomes), function(x) {
		n = chromosomes[x]
		if (n > 0) {
			gff.chr <- gff.gr[gff.gr@seqnames %in% x]
			sample(gff.chr, n, replace = F)
		}
	})
}


genomic.clustering <- function(list.genes, gff.gr, random.samplings) {
	gr.subset <- gff.gr[elementMetadata(gff.gr)[ ,'gene_name'] %in% list.genes]
	distance.ranges <- gr.dist(gr.subset) %>% .[upper.tri(., diag = F)]
	distance.ranges <- distance.ranges[!is.na(distance.ranges)]
	
	chromosomes <- table(gr.subset@seqnames)

	set.seed(123)
	sampled.distances <- list()
	for (index in 1:random.samplings) {
		sampled.gr <- sample.genes(chromosomes, gff.gr)	%>% plyranges::bind_ranges()
		distance.sampled <- gr.dist(sampled.gr) %>% .[upper.tri(., diag = F)]
		distance.sampled <- distance.sampled[!is.na(distance.sampled)]
		list.name <- paste("random", index, sep = "")
		sampled.distances[[list.name]] <- distance.sampled
		message(glue("computed distances for random sampling # {index}"))
	}
	return(list(actual.dist = distance.ranges, random.dist = sampled.distances, chr.distribution = chromosomes))
}

gff = "~/Downloads/gencode.v43.annotation.gff3"
gff.gr = rtracklayer::import(gff) # creates a GRanges object
gff.gr <- gff.gr[elementMetadata(gff.gr)[ ,'type'] == "gene"]


cellTypes <- names(signature.list)
direction.signatures <- "positive"
distance.genes.cellTypes <- list()
for (cellType in cellTypes) {
	geneset <- signature.list[[cellType]]
	genes.dist.list <- list()
	for (dims in names(geneset)) {
		l1 <- geneset[[dims]][[direction.signatures]]
		if (length(l1) > 2) {
			#geneset.new[[dims]] <- l1
			distances.genes <- genomic.clustering(l1, gff.gr, random.samplings = 10)
			
			####in cases were none of the signature genes are present on same chromosome
			if (length(distances.genes[['actual.dist']] > 1)) {
				actual.dist <- distances.genes[['actual.dist']] %>% data.frame("iteration" = "Actual", distance = .)
				random.dist <- distances.genes[['random.dist']] %>% lapply(., function(x) data.frame(distance = x))
				random.dist <- bind_rows(random.dist, .id = "iteration")
			
				genes.dist.list[[dims]] <- rbind(actual.dist, random.dist)
			}
			
		}
	}
	genes.dist.list <- bind_rows(genes.dist.list, .id = "Dim")
	distance.genes.cellTypes[[cellType]] <- genes.dist.list
	
	message("*************")
	message("*************")
	message(glue("done for {cellType}"))
}

distance.genes.cellTypes <- bind_rows(distance.genes.cellTypes, .id = "cellType")
distance.genes.cellTypes$type <- ifelse(distance.genes.cellTypes$iteration == "Actual", "Actual", "Random")

write.table(distance.genes.cellTypes, file = glue("genomic_distances_ICs.{direction.signatures}Genes.txt"), sep = "\t", quote = F, row.names = F)




######make plots

direction.signatures <- "negative"
data.plot <- read.table(glue("genomic_distances_ICs.{direction.signatures}Genes.txt"), header = T)
cellTypes <- unique(data.plot$cellType)
list.plots <- list()
for (cellType in cellTypes) {
	data.plot <- distance.genes.cellTypes[distance.genes.cellTypes$cellType == cellType, ]
	
	graphic <- ggplot(data = data.plot, aes(x = type, y = distance)) + geom_boxplot() + 
			xlab("") + ylab("#base pairs") + ggtitle(cellType) + facet_wrap(~Dim, nrow = 1, strip.position = "right") +
			stat_compare_means(method = "wilcox", label = "p.format", label.y = 2.0e+8, label.x.npc = "left", size = 3, col = "darkred")
			
	#list.plots[[cellType]] <- graphic + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 25, vjust = 1, hjust = 1))
	list.plots[[cellType]] <- graphic + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), plot.margin = unit(c(0, 0, -0.5, 0), "cm"))
}


file.name <- glue("genomicDistances.ICA.signatureGenes-{direction.signatures}.{tolower(cohort)}.IDH{idh}--filtered-noNorm-optimalRank.svg")

svglite(file.name, height = 8)
do.call("grid.arrange", c(list.plots, ncol = 1)) 
dev.off()

