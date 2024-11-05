library(dplyr); library(magrittr);
library(glue);

source('../../melanomaSplicing/generalFunctions.R')

get.genes <- function(index, S) {
	x <- S[ ,index]
	pos <- x[x > mean(x) + list.cutoffs[index]] %>% names()
	neg <- x[x < mean(x) - list.cutoffs[index]] %>% names()
	genes <- list("positive" = pos, "negative" = neg)
	return(genes)
}

genesets.ica <- function(ica.res, cutoff = 2.5) {
	S <- ica.res$S
	stdev <- apply(S, 2, sd)
	
	list.cutoffs <- cutoff * stdev
	
	gene.sets <- sapply(1:length(list.cutoffs), function(index) get.genes(index, S), simplify = F)
	return(gene.sets)
}
#####select only the

sampleType = "mut"
cohort <- "TCGA"
dims <- 10
load(glue("ICA-JADE-{dims}.cellTypes.{tolower(cohort)}.IDH{sampleType}-filtered-noNorm-optimalRank.Rda"))


signature.list <- list()
for(cellType in names(ica.res.list)) {
	ica.res <- ica.res.list[[cellType]]
	gene.list <- genesets.ica(ica.res, cutoff = 2.5) %>% set_names(., paste("Dim", 1:length(.), sep = "."))
	signature.list[[cellType]] <- gene.list
}
save(signature.list, file = glue("ICA-JADE-{dims}.signatureGenes.{tolower(cohort)}.IDH{sampleType}--filtered-noNorm-optimalRank.Rda"))
