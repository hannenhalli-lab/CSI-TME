library(MineICA); library(dplyr)
library(magrittr); library(glue)
library(doParallel); library(doSNOW)

#source('clusterICA.modified.R')

preProcessInput <- function(inputFile, removeNA = T, removeSD = T, removeMedian = T, normalize = F, logT = F, zscore = F, positive = F) {
	if (grep("[A-Z]", inputFile[ ,1][1])) {
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
source('functions.CSI-TME.R')
cohort.discovery <- "tcga"
#baseDir <- "~/DataAndAnalysis/cellTypeDeconvolution/PCA"
baseDir <- "../data/"
inputPath <- glue("{baseDir}/{cohort.discovery}_{idh}_filtered/")
cellTypes <- list.files(inputPath) %>% gsub(".txt", "", .)
outputPath <- "./"

dims = "10"; max.dim <- 10; method = "JADE";

output.name <- glue("ICA-{method}-{dims}.cellTypes.{cohort.discovery}.IDH{idh}-filtered-noNorm-optimalRank.Rda")




###progress monitor solution from  https://stackoverflow.com/questions/5423760/how-do-you-create-a-progress-bar-when-using-the-foreach-function-in-r


dims.cellTypes <- c(Bcell = 42, Endothelial = 40, Malignant = 59, Myeloid = 47, Oligos = 52, Stromal = 44, Tcell = 45)
ica.res.list <- list()
for (cellType in cellTypes) {
	inputFile <- read.table(glue("{inputPath}/{cellType}.txt"), sep = "\t", header = T)
	expression.data <- preProcessInput(inputFile, removeNA = T, removeSD = T, removeMedian = T, logT = T, normalize = F, zscore = T, positive = F)
	if (dims == "MSTD") {
		max.dim <- dims.cellTypes[[cellType]]
	}
	if (method == "JADE") {
		res <- runICA(method = "JADE", X = expression.data, nbComp = max.dim, alg.type = "parallel", maxit = 20000, tol = 10^-6)
			
	} else {
		cl <- makeCluster(future::availableCores() - 4)
		registerDoSNOW(cl)
		res <- clusterFastICARuns.modified(X = expression.data, nbComp = max.dim, alg.type = "deflation", nbIt = 50, funClus = "hclust", method = "average")
	
	}
	ica.res.list[[cellType]] <- res
	print(glue("done for {cellType}"))
}
stopCluster(cl);

save(ica.res.list, file = output.name)


