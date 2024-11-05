#library(clusterProfiler)
library(dplyr); library(magrittr);
library(gtools); library(preprocessCore); 
library(reshape2); library(ggplot2);
library(glue); library(TCGAbiolinks); library(purrr);
library(survival); library(doParallel)


#####select only the

sampleType = "mut"
source('../functionsForInteractionAnalysis.R')
cohort <- "TCGA"; cohort.names <- "tcga"
method <- "JADE"; dims = "10" ####for dim 25, pvalue of 0.006 corresponds to the FDR of 0.01024762, for dim 10, pvalue of 0.007 is FDR of 1 perecnt

load(glue("ICA-{method}-{dims}.cellTypes.{tolower(cohort)}.IDH{sampleType}-filtered-noNorm-optimalRank.Rda"))

binMapList <- componentBinningFunction(pcaList = ica.res.list, factorization = "ICA", sample.names = NULL)

#binMapList <- componentBinningFunction(pcaList = ica.res.list, factorization = "NMF", sample.names = sample.names[-1])


if (cohort == "TCGA") {
	cancerType <- c("GBM", "LGG")
	clinical <- read.table("/Users/singha30/DataAndAnalysis/SPAGE/RB1-Project/AveragedSamples/AllAvialableClinicalDataMod.txt", sep = "\t", header = T)
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

cellTypeInteractionList <- list()
progressIndex = 0
covariateType <- "technical" ## or "biological"
#covariateType <- "biological" ## or "technical"


interactionFile <- read.table(glue("ICA-JADE-{dims}.allInteractions-IDHmut-TCGA-technical.txt"), sep = "\t", header = T)
interactionFile <- na.omit(interactionFile)
interactionFile <- interactionFile[interactionFile$FDR < 0.25, ]
#nteractionFile <- interactionFile[interactionFile$pvalue < 0.007, ]

iteration.stats.list <- list()
for (row.index in 1:nrow(interactionFile)) {
	int.line <- interactionFile[row.index, ]
	sign.int <- sign(int.line$HR)
	#if (int.line$FDR < 0.30) {
		cellTypes <- strsplit(int.line$cellTypes, split = ":") %>% unlist()
		components <- strsplit(int.line$components, split = ":") %>% unlist()
		cellType1 <- cellTypes[1]; component1 = components[1]
		cellType2 <- cellTypes[2]; component2 = components[2]
		int.bin <- int.line[ ,'bin']
		
		score.1 <- ica.res.list[[cellType1]]$A %>% set_colnames(paste("Dim", 1:dims, sep = "."))
		score.2 <- ica.res.list[[cellType2]]$A %>% set_colnames(paste("Dim", 1:dims, sep = "."))
		
		comp1.values <- score.1[ ,component1]; comp1.bins <- binMapList[[cellType1]][ ,component1]
		comp2.values <- score.2[ ,component2]; comp2.bins <- binMapList[[cellType2]][ ,component2]
		
		if (all(rownames(comp1.values) == rownames(comp2.values))) {
			componentvalues <- data.frame(samples = names(comp1.values), comp1 = comp1.values, comp2 = comp2.values) %>% set_rownames(NULL)
			binValues <- data.frame(samples = names(comp1.bins), bin1 = comp1.bins, bin2 = comp2.bins) %>% set_rownames(NULL)
			
			componentvalues$samples = gsub("\\.", "-", componentvalues$samples) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
			binValues$samples = gsub("\\.", "-", binValues$samples) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
			
			if (int.bin == "Bin.1") {
				bins <- ifelse(binValues$bin1 == 0 & binValues$bin2 == 0, 1, 0)
			}
			if (int.bin == "Bin.3") {
				bins <- ifelse(binValues$bin1 == 2 & binValues$bin2 == 0, 1, 0)
			}
			if (int.bin == "Bin.9") {
				bins <- ifelse(binValues$bin1 == 2 & binValues$bin2 == 2, 1, 0)
			}
			
			variables <- cbind(componentvalues, bins)
			clinicalData <- merge(clinical, variables, by = "samples")
			clinicalData <- clinicalData[!duplicated(clinicalData$samples), ]
			
			iter.values <- list()
			for (iter in 1:11) {
				if (iter == 1) {
					indexes <- sample(1:nrow(clinicalData), nrow(clinicalData), replace = F)
					dt1 <- clinicalData[indexes, ]
				} else {
					size <- floor(nrow(clinicalData) * 80 / 100)
					indexes <- sample(1:nrow(clinicalData), size, replace = F)
					dt1 <- clinicalData[indexes, ]
				}
				tryCatch({
					cox.out = coxph(Surv(time,status) ~ bins + comp1 + comp2 + age + sex, data = dt1)
					coefficients = summary(cox.out)$coefficients
					hr <- coefficients['bins', 1]
					pvalue <- coefficients['bins', 5]
					values <- data.frame(hr = hr, pvalue = pvalue)
					iter.values[[iter]] <- values
			    },error=function(x){})
			}
			iter.values <- bind_rows(iter.values, .id = "iteration")
			iter.stats <- sum(sign(iter.values$hr) == sign.int & iter.values$pvalue < 0.05)
			iteration.stats.list[[row.index]] <- iter.stats
		}
		#}
	message(glue("done for {row.index} interactions"))
}

interactionFile$Reproduciblity <- unlist(iteration.stats.list)
#filtered.interactions <- interactionFile[interactionFile$Reproduciblity >= 8, ]
write.table(interactionFile, file = glue("ICA-{method}-{dims}.FDR0.25.robustInteractions-IDH{sampleType}-{cohort}-{covariateType}.txt"), sep = "\t", quote = F, row.names = F)


