#library(clusterProfiler)
library(dplyr); library(magrittr);
library(preprocessCore); 
library(glue); library(TCGAbiolinks); library(purrr);
library(survival); library(doParallel);


preProcessInput <- function(inputFile, removeNA = T, removeSD = T, removeMedian = T, normalize = F, logT = F, zscore = F, positive = F) {
	if (length(grep("[A-Z]", inputFile[ ,1][1])) == 1) {
		inputFile <- inputFile[ ,-1] %>% set_rownames(inputFile[ ,1])
	}
	if (removeNA) {
		index1 <- apply(inputFile, 1, function(x) sum(is.na(x))) == ncol(inputFile)
		subData <- inputFile[!index1, ]
	}else{
		subData <- inputFile
	}
	if (removeSD) {
		index2 <- apply(subData, 1, function(x) sd(x, na.rm = T)) == 0
		subData <- subData[!index2, ]
	}
	if (removeMedian) {
	    index3 <- apply(subData, 1, function(x) median(x, na.rm = T)) < 1
	    subData <- subData[!index3, ]
	}
    if (logT) {
        subData <- log2(subData + 0.001)
    }
	if (normalize) {
		subData <- subData %>% data.matrix %>% normalize.quantiles %>% 
			set_rownames(rownames(subData)) %>% set_colnames(colnames(subData))
	}
	if (zscore == T) {
		subData <- subData %>% apply(., 1, scale) %>% t() %>% set_colnames(colnames(subData))
		if (positive == T) {
			subData <- nneg(subData, "min")
		}
	}
	#subData <- subData %>% t()
	return(subData)
}



idh = "mut"
source('../functionsForInteractionAnalysis.R')
cohort.discovery <- "TCGA"
#cohort.validation <- "immunotherapy"
#baseDir <- "~/DataAndAnalysis/cellTypeDeconvolution/PCA"

cohort.validation <- "IDHtreatment"
baseDir <- "IDHinhibitor"


inputPath <- glue("{baseDir}/{tolower(cohort.validation)}_{idh}_filtered/")
outputPath <- "./"
dims = "10"

load(glue("ICA-JADE-{dims}.cellTypes.{tolower(cohort.discovery)}.IDH{idh}-filtered-noNorm-optimalRank.Rda"))
discovery.pca <- ica.res.list
output.name <-glue("ICA-JADE-{dims}.cellTypes.{tolower(cohort.validation)}.IDH{idh}-filtered-noNorm-optimalRank.Rda")
if (cohort.validation == "immunotherapy") {
	inputPath <- glue("{baseDir}/{tolower(cohort.validation)}_glioma_filtered/")
	output.name <- glue("ICA-JADE-{dims}.cellTypes.{tolower(cohort.validation)}.IDHglioma-filtered-noNorm-optimalRank.Rda")
}

ica.res.list <- list(); remove.median = T ##do not filter
for (cellType in names(discovery.pca)) {
	pca1 <- discovery.pca[[cellType]]
	
	inputFile <- read.table(glue("{inputPath}/{cellType}.txt"), sep = "\t", header = T)
	expression.data <- preProcessInput(inputFile, removeNA = T, removeSD = T, removeMedian = remove.median, logT = T, normalize = F, zscore = F, positive = F)
	
	discovery.rotation  <- pca1$S
	
	commonGenes <- rownames(expression.data)[rownames(expression.data) %in% rownames(discovery.rotation)]
	
	discovery.rotation <- discovery.rotation[commonGenes, ]
	expression.data <- expression.data[commonGenes, ]
	
	scores.projection <- scale(t(expression.data), center = T, scale = T) %*% discovery.rotation
	scores.projection <- scores.projection / 1000 ####scale down
	ica <- list()
	ica$A <- scores.projection
	ica.res.list[[cellType]] <- ica
	
	message(glue("done for {cellType}"))
}
save(ica.res.list, file = output.name)


######CROSS cell types#####

output.name <-glue("ICA-JADE-{dims}.crossCellTypes.{tolower(cohort.validation)}.IDH{idh}-filtered-noNorm-optimalRank.Rda")

cellTypeCombinations <- combn(names(ica.res.list), 2)


expression.files <- c("Bcell.txt","Endothelial.txt","Malignant.txt","Myeloid.txt","Oligos.txt","Stromal.txt","Tcell.txt")
expression.list <- list()
for (expression.file in expression.files) {
	cellType <- gsub(".txt", "", expression.file)
	componentMatrix1 <- read.table(glue("{inputPath}/{expression.file}"), header = T)
	componentMatrix1 <- preProcessInput(componentMatrix1, removeNA = T, removeSD = T, removeMedian = T, logT = T, normalize = F, zscore = F, positive = F)
	expression.list[[cellType]] <- 	componentMatrix1
	message(glue("read expression file for {expression.file}"))
}


ica.res.list <- list()
for (cellTypeIndex in 1:ncol(cellTypeCombinations)) {
	cellType1 <- cellTypeCombinations[1, cellTypeIndex]
	cellType2 <- cellTypeCombinations[2, cellTypeIndex]

	pca <- discovery.pca[[cellType1]] ##discovery
	expression.data <- expression.list[[cellType2]] ##validation
	discovery.rotation  <- pca$S
		
	commonGenes <- rownames(expression.data)[rownames(expression.data) %in% rownames(discovery.rotation)]
	
	discovery.rotation <- discovery.rotation[commonGenes, ]
	expression.data <- expression.data[commonGenes, ]
	
	
	scores.projection <- scale(t(expression.data), center = T, scale = T) %*% discovery.rotation
	scores.projection <- scores.projection / 1000 ####scale down
	ica <- list()
	ica$A <- scores.projection
	
	cellType.pair <- paste(cellType1, ":", cellType2, sep = "")
	ica.res.list[[cellType.pair]] <- ica
	
	message(glue("done for {cellTypeIndex} celltype combinations"))
}

save(ica.res.list, file = output.name)

