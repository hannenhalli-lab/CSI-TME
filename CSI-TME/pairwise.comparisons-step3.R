#library(clusterProfiler)
library(dplyr); library(magrittr);
#library(gtools); library(preprocessCore); 
library(reshape2); library(ggplot2);
library(glue); library(TCGAbiolinks);
library(survival); library(doParallel)


#####select only the

estimateFDR <- F; niter = 10 ####only used when bootstrapping 

sampleType = "mut"
source('../functionsForInteractionAnalysis.R')
cohort <- "TCGA"; cohort.names <- "tcga"
method <- "JADE"; dims = "10"

load(glue("ICA-{method}-{dims}.cellTypes.{tolower(cohort)}.IDH{sampleType}-filtered-noNorm-optimalRank.Rda"))

binMapList <- componentBinningFunction(pcaList = ica.res.list, factorization = "ICA", sample.names = NULL)

#binMapList <- componentBinningFunction(pcaList = ica.res.list, factorization = "NMF", sample.names = sample.names[-1])


if (cohort == "TCGA") {
	cancerType <- c("GBM", "LGG")
	clinical <- read.table("AllAvialableClinicalDataMod.txt", sep = "\t", header = T)
	clinical <- clinical[clinical$type %in% cancerType, ]
	
	codel <- read.csv("../All_metadata-CCGA.csv", header = T) %>% set_names(c('samples' ,'project', 'time', 'status', 'IDH', 'grade', 'age', 'subtype', 'codel'))
	codel <- codel[grep("TCGA", codel$samples), c('IDH', 'samples', 'codel', 'grade')]
	codel$samples <- gsub("\\.", "-", codel$samples) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
	clinical <- merge(clinical, codel, by = "samples")
} else {
	if (grep("CGGA", cohort)) {
		clinical <- read.csv("../All_metadata-CCGA.csv", header = T) %>% set_names(c('samples' ,'project', 'time', 'status', 'IDH', 'grade', 'age', 'subtype', 'codel'))
		gender <- read.table("../PCA/cgga.combinedMetadata.txt", sep = "\t", header = T) %>% .[ ,c('samples', 'Gender')]
		clinical <- merge(clinical, gender, by = "samples")
		names(clinical)[names(clinical) == 'Gender'] <- "sex"
	}
}



cellTypeCombinations <- combn(names(binMapList), 2)

####parallel solution adapted from https://stackoverflow.com/questions/38318139/run-a-for-loop-in-parallel-in-r#######
#symmetricBins <- c(4,7,8)

cellTypeInteractionList <- list(); iteraction.iterations <- list()
progressIndex = 0
covariateType <- "technical" ## or "biological"
#covariateType <- "biological" ## or "technical"

set.seed(123)
for (iteration in 1:niter) {
	if (estimateFDR == TRUE) {
		clinical$samples <-  sample(clinical$samples, replace = F)
	}
	
	for (cellTypeIndex in 1:ncol(cellTypeCombinations)) {
		progressIndex <- progressIndex + 1
	
		cellType1 <- cellTypeCombinations[1, cellTypeIndex]
		cellType2 <- cellTypeCombinations[2, cellTypeIndex]
	
		componentMatix1 <- ica.res.list[[cellType1]]$A
		componentMatix2 <- ica.res.list[[cellType2]]$A
	
		ncomp1 <- ncol(componentMatix1)
		ncomp2 <- ncol(componentMatix2)
	
	
		colnames(componentMatix1) <- paste("Dim.", 1:ncomp1, sep = "")
		colnames(componentMatix2) <- paste("Dim.", 1:ncomp2, sep = "")
	
	
		componentInteractionList <- NA
	
		binMap1 <- binMapList[[cellType1]] #%>% set_rownames(samplenames1)
		binMap2 <- binMapList[[cellType2]] #%>% set_rownames(samplenames1)
	
		if (is.null(ncol(binMap1)) | is.null(ncol(binMap2))) {
			cellTypeInteractionList[[cellTypeIndex]] <- componentInteractionList
			next
		}else {
		
		
			if (all(rownames(binMap1) == rownames(binMap2))) {
			
				if (cohort == "TCGA") {
					allSampleNames <- gsub("\\.", "-", rownames(binMap1))
					rownames(binMap1) <- allSampleNames
					rownames(binMap2) <- allSampleNames
					rownames(componentMatix1) <- allSampleNames
					rownames(componentMatix2) <- allSampleNames
				
					desiredSamples <- TCGAquery_SampleTypes(allSampleNames, typesample = c("TP", "TR", "TM"))
				} else {
					if (grep("CGGA", cohort)) {
						allSampleNames <- rownames(binMap1)
						desiredSamples <- grep("CGGA", allSampleNames, value = T)
					}
				}
			
			
			
			}
	
		
			binMap1 <- binMap1[rownames(binMap1) %in% desiredSamples, 1:ncomp1]
			binMap2 <- binMap2[rownames(binMap2) %in% desiredSamples, 1:ncomp2] 
		
			componentMatix1 <- componentMatix1[rownames(componentMatix1) %in% desiredSamples, 1:ncomp1]
			componentMatix2 <- componentMatix2[rownames(componentMatix2) %in% desiredSamples, 1:ncomp2] 
			
			#componentMatix1 <- componentMatix1[ ,!colnames(componentMatix1) %in% badComponents1]
			#componentMatix2 <- componentMatix2[ ,!colnames(componentMatix2) %in% badComponents2]
			
			#binMap1 <- binMap1[ ,!colnames(binMap1) %in% badComponents1]
			#binMap2 <- binMap2[ ,!colnames(binMap2) %in% badComponents2]
		
		
			if (cohort == "TCGA") { 
				rownames(binMap1) <- rownames(binMap1) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
				rownames(binMap2) <- rownames(binMap2) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
				rownames(componentMatix1) <- rownames(componentMatix1) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
				rownames(componentMatix2) <- rownames(componentMatix2) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
			}
		
				#bin.map <- cbind(binMap1[ ,component1], binMap2[ ,component2])
		
		
		
			lv.1 <- colnames(componentMatix1); lv.2 <- colnames(componentMatix2)
		
			componentCombinations <- expand.grid(lv.1, lv.2)
			sizeCount <- nrow(componentCombinations)
			cores=detectCores()
			cl <- makeCluster(cores[1] - 6) #not to overload your computer
			#cl <- makeCluster(9)
			registerDoParallel(cl)
			#for (componentIndex in 1:ncol(componentCombinations)) {
			componentInteractionList <- NULL
			componentInteractionList <- foreach(componentIndex=1:sizeCount, .packages=c('survival', 'dplyr', 'magrittr')) %dopar% 
					{
						component1 <- componentCombinations[componentIndex, 1]
						component2 <- componentCombinations[componentIndex, 2]
					
					
						bin.map <- cbind(binMap1[ ,component1], binMap2[ ,component2])
					    #bins = bin.map[ ,1] * 3 + bin.map[ ,2] + 1 #####this one is used in spageFinder
						bins = bin.map[ ,2] * 3 + bin.map[ ,1] + 1 ####I switched the indexes for calculation to call bins in order of comp1:comp2
					  	bins.binary = sapply(1:9, function(tt) ifelse(bins == tt, 1, 0))
					    colnames(bins.binary) = sapply(1:9,function(v) paste('bin',v,sep=''))
					
						comp1 <- componentMatix1[ ,component1] #%>% set_names(samplenames)
						comp2 <- componentMatix2[ ,component2] #%>% set_names(samplenames)
					
						componentvalues <- data.frame(samples = names(comp1), comp1 = comp1, comp2 = comp2)
					
						clinicalData <- merge(clinical, componentvalues, by = 'samples')
					  	#bin.binary = bins.binary[,bin]
						#survResults <- foreach(binIndex=1:sizeCount, .packages=c('survival', 'dplyr')) %dopar%
						tempList <- list()
						#nameInfo <- list()
						#desiredComponents <- c(1, 3, 9)
						desiredBins <- c(1,3,7,9)
						for (binIndex in desiredBins) {
							tempMatrix = perform.cox(binaryBins = bins.binary, index = binIndex, clinicalData = clinicalData, type = sampleType, covariate = covariateType) #calling a function
							tempList[[binIndex]] <- tempMatrix
							#nameInfo[[binIndex]] <- c(component1, component2, binIndex)
						}
						tempList #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
							#}
					}
		
			#}
			#stop cluster
			stopCluster(cl)
			componentNames <- glue("{componentCombinations[, 1]}:{componentCombinations[, 2]}")
			names(componentInteractionList) <- componentNames
			cellTypeInteractionList[[cellTypeIndex]] <- componentInteractionList
		}
		print(glue("done for {cellTypeIndex} cellType combinations"))
	}
	names(cellTypeInteractionList) <- glue("{cellTypeCombinations[1, ]}:{cellTypeCombinations[2, ]}")
	
	if (estimateFDR == FALSE) {
		break
	} else {
		iteraction.iterations[[iteration]] <- cellTypeInteractionList
	}
}

binNames <- paste("Bin", c(1,3,7,9), sep = ".")
#componentLength <- ncol(componentCombinations)
HRparameters <- list(
	#componentLength = componentLength, 
	selection = 1, 
	#componentNames = componentNames, 
	binNames = binNames
	)
	
Pvalueparameters <- list(
	#componentLength = componentLength, 
	selection = 2, 
	#componentNames = componentNames, 
	binNames = binNames
	)
#cellTypeInteractionDf <- outputParser(inputList = cellTypeInteractionList[1], parameterList = parameterList)



if (estimateFDR == F) {
	cellTypeInteractionList <- cellTypeInteractionList[!is.na(cellTypeInteractionList)]
	subString <- "CIBERSORTxHiRes_Job1_|Window\\w+\\.txt|_"
	hazardRatiosDf <- NULL
	pvalueDf <- NULL
	for (listIndex in 1:length(cellTypeInteractionList)) {
		cellTypes <- names(cellTypeInteractionList)[listIndex] %>% gsub(glue(subString), "", .)
		inputList <- cellTypeInteractionList[listIndex]
		tempList <- outputParser(inputList = inputList, parameterList = HRparameters) %>% 
				mutate( cellTypes = cellTypes, components = rownames(.), .before = "Bin.1")
			
		hazardRatiosDf <- rbind(hazardRatiosDf, tempList)
	
		tempList <- outputParser(inputList = inputList, parameterList = Pvalueparameters) %>% 
				mutate( cellTypes = cellTypes, components = rownames(.), .before = "Bin.1")
			
		pvalueDf <- rbind(pvalueDf, tempList)
	}

	significantInteractions <- getSignificantInteractions(hazardRatiosDf = hazardRatiosDf, pvalueDf = pvalueDf, method = "fdr", threshold = 1, factorList = "cellTypes", keep.bins = c("Bin.1", "Bin.3", "Bin.7", "Bin.9"), symmetricBinCorrection = T)
	write.table(significantInteractions, file = glue("ICA-{method}-{dims}.allInteractions-IDH{sampleType}-{cohort}-{covariateType}.txt"), sep = "\t", quote = F, row.names = F)
	
}

#binNames <- paste("Bin", c(1,3,7,9), sep = ".")


significant.list <- list()
if (estimateFDR == T) {
	for (iteration in 1:niter) {
		cellTypeInteractionList <- iteraction.iterations[[iteration]]
		cellTypeInteractionList <- cellTypeInteractionList[!is.na(cellTypeInteractionList)]
		subString <- "CIBERSORTxHiRes_Job1_|Window\\w+\\.txt|_"
		hazardRatiosDf <- NULL
		pvalueDf <- NULL
		for (listIndex in 1:length(cellTypeInteractionList)) {
			cellTypes <- names(cellTypeInteractionList)[listIndex] %>% gsub(glue(subString), "", .)
			inputList <- cellTypeInteractionList[listIndex]
			tempList <- outputParser(inputList = inputList, parameterList = HRparameters) %>% 
					mutate( cellTypes = cellTypes, components = rownames(.), .before = "Bin.1")
			
			hazardRatiosDf <- rbind(hazardRatiosDf, tempList)
	
			tempList <- outputParser(inputList = inputList, parameterList = Pvalueparameters) %>% 
					mutate( cellTypes = cellTypes, components = rownames(.), .before = "Bin.1")
			
			pvalueDf <- rbind(pvalueDf, tempList)
		}

		significantInteractions <- getSignificantInteractions(hazardRatiosDf = hazardRatiosDf, pvalueDf = pvalueDf, method = "fdr", threshold = 1, factorList = "cellTypes", keep.bins = c("Bin.1", "Bin.3", "Bin.7", "Bin.9"), symmetricBinCorrection = T)	
		significant.list[[iteration]] <- significantInteractions
	}
	#
}


significantInteractions <- bind_rows(significant.list, .id = "iteration")
write.table(significantInteractions, file = glue("ICA-{method}-{dims}.randomInteractions-IDH{sampleType}-{cohort}-{covariateType}.txt"), sep = "\t", quote = F, row.names = F)
