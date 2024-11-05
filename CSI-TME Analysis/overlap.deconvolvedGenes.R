library(magrittr); library(glue); library(dplyr)
library(preprocessCore)
orgDb <- "org.Hs.eg.db"

source('../../melanomaSplicing/generalFunctions.R')
source('../functionsForInteractionAnalysis.R')


expression.path.tcga <- "tcga_mut_filtered/"
expression.path.cgga <- "cgga693_mut_filtered/"

files <- list.files(expression.path.tcga)

get.topGenes <- function(gene.vector, n = NULL, p = NULL) {
	if (is.null(n)) {
		size.x = length(gene.vector); 
		top.genes.n = floor((size.x)* p) / 100; 
		gene.list <- gene.vector[1:top.genes.n] %>% names()
	}
	if (is.null(p)) {
		gene.list <- gene.vector[1:top.genes.n] %>% names()
	}
	return(gene.list)
}

get.jaccard <- function(a, b) {
	###from https://www.r-bloggers.com/2021/11/how-to-calculate-jaccard-similarity-in-r-2/#####
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}

jaccard.test.wrapper <- function(a, b, output) {
	set.seed(1234)
	total <- union(a, b)
	vec1 <- ifelse(total %in% a, 1, 0)
	vec2 <- ifelse(total %in% b, 1, 0)
	j.test <- jaccard.test(vec1, vec2, method = "bootstrap")
	pvalue <- j.test$pvalue
	if (output == "pvalue") {
		return(pvalue)
	}
	if (output == "object") {
		return(j.test)
	}
}

tcga.expression.list <- list()
for (file.name in files) {
	gene.expression <- read.table(glue("{expression.path.tcga}/{file.name}"), header = T, row.names = 1)
	gene.expression <- log2(gene.expression + 0.001)
	gene.expression <- normalize.quantiles(data.matrix(gene.expression)) %>% 
		set_rownames(rownames(gene.expression)) %>% set_colnames(colnames(gene.expression))
	
	cellType <- gsub("\\.txt", "", file.name)
	tcga.expression.list[[cellType]] <- gene.expression
	message(glue("read data for {file.name}"))
}

cgga.expression.list <- list()
for (file.name in files) {
	gene.expression <- read.table(glue("{expression.path.cgga}/{file.name}"), header = T, row.names = 1)
	gene.expression <- log2(gene.expression + 0.001)
	gene.expression <- normalize.quantiles(data.matrix(gene.expression)) %>% 
		set_rownames(rownames(gene.expression)) %>% set_colnames(colnames(gene.expression))
	
	cellType <- gsub("\\.txt", "", file.name)
	cgga.expression.list[[cellType]] <- gene.expression
	message(glue("read data for {file.name}"))
}



tcga.median.list <- list(); cgga.median.list <- list()
for (cellType in names(tcga.expression.list)) {
	tcga.median <- apply(tcga.expression.list[[cellType]], 1, median, na.rm = T)
	tcga.median <- tcga.median[order(tcga.median, decreasing = T)]
	
	cgga.median <- apply(cgga.expression.list[[cellType]], 1, median, na.rm = T)
	cgga.median <- cgga.median[order(cgga.median, decreasing = T)]
	
	tcga.median.list[[cellType]] <- tcga.median
	cgga.median.list[[cellType]] <- cgga.median
}

tcga.gene.list <- lapply(tcga.median.list, get.topGenes, p = 10)
cgga.gene.list <- lapply(cgga.median.list, get.topGenes, p = 10)

gene.frequency <- unlist(tcga.gene.list) %>% table
gene.specific <- gene.frequency[gene.frequency < 3] %>% names()
gene.ubiqutous <- gene.frequency[gene.frequency >= 3] %>% names()
tcga.gene.list <- lapply(tcga.gene.list, function(x) x[x %in% gene.specific])

gene.frequency <- unlist(cgga.gene.list) %>% table
gene.specific <- gene.frequency[gene.frequency < 3] %>% names()
gene.ubiqutous <- gene.frequency[gene.frequency >= 3] %>% names()
cgga.gene.list <- lapply(cgga.gene.list, function(x) x[x %in% gene.specific])

JI.list <- list(); test = "fisher"

for (cellType.tcga in names(tcga.gene.list)) {
	tcga.list <- tcga.gene.list[[cellType.tcga]]
	tcga.universe <- rownames(tcga.expression.list[[cellType.tcga]])
	ji.list <- list()
	for (cellType.cgga in names(cgga.expression.list)) {
		cgga.list <- cgga.gene.list[[cellType.cgga]]
		
		cgga.universe <- rownames(cgga.expression.list[[cellType.cgga]])
		#universe <- tcga.universe[tcga.universe %in% cgga.universe]
		universe <- intersect(tcga.universe, cgga.universe)
		if (test == "fisher") {
			fisher.res <- fisherTestFunction(tcga.list, cgga.list, universe = universe)
			ji.list[[cellType.cgga]] <- data.frame(CGGA = cellType.cgga, oddsRatio = fisher.res[4], pvalue = fisher.res[5])
		} else {
			JI <- get.jaccard(tcga.list, cgga.list);
			#pvalue <- jaccard.test.wrapper(a = tcga.list, b = cgga.list, output = "pvalue") ##package doesnt work anymore
			pvalue <- NA
			ji.list[[cellType.cgga]] <- data.frame(CGGA = cellType.cgga, JI = JI, pvalue = pvalue)
		}
	}
	JI.list[[cellType.tcga]] <- do.call("rbind", ji.list)
}

JI.data <- bind_rows(JI.list, .id = "TCGA")
#data.plot <- reshape2::dcast(data = JI.data,formula = TCGA ~ CGGA, value.var = "JI")
#data.plot <- data.plot[ ,-1] %>% set_rownames(data.plot[ ,1])
#pheatmap::pheatmap(data.plot, scale = 'row')

data.plot <- JI.data
data.plot$type <- ifelse(JI.data$TCGA ==  JI.data$CGGA, "same cell type", "across cell types")
data.plot$type <- factor(data.plot$type, levels = c("same cell type", "across cell types"))


write.table(data.plot, file = "overlap.deconvovledTopGenes.cellTypes.txt", sep = "\t", quote = F, row.names = F)

if (test == "JI") {
	output.file <- glue("plots/boxplot.jaccard.topGenes.TCGA.CGGA.svg")
	svglite::svglite(file = output.file, height = 2.5, width = 2.5)
	ggplot(data = data.plot, aes(x = type, y = JI)) + geom_boxplot() + theme_bw() + ggpubr::stat_compare_means(label = "p.format", label.x = 1.5)
	dev.off()
}

if (test == "fisher") {
	data.plot$FDR <- p.adjust(data.plot$pvalue, method = "fdr")
	output.file <- glue("plots/boxplot.overlap.topGenes.TCGA.CGGA.svg")
	svglite::svglite(file = output.file, height = 2.5, width = 2.5)
	ggplot(data = data.plot, aes(x = type, y = oddsRatio)) + geom_boxplot() + theme_bw() + ggpubr::stat_compare_means(label = "p.format")
	dev.off()
}


data.plot <- data.plot[data.plot$type == "same cell type", ]
data.plot$signi <- -log10(data.plot$FDR)

output.file <- glue("plots/dotplot.overlap.topGenes.TCGA.CGGA-SameCellType.svg")
svglite::svglite(file = output.file, height = 3.5, width = 2.5)
ggplot(data = data.plot, aes(x = "overlap", y = TCGA)) + geom_point(aes(size=oddsRatio)) + scale_size_continuous(range = c(2,5), breaks = c(2,8,14)) + theme_bw()
dev.off()





