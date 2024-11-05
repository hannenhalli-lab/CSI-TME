#library(clusterProfiler)
library(dplyr); library(magrittr);
library(gtools); library(preprocessCore); 
library(reshape2); library(ggplot2);
library(glue); library(TCGAbiolinks); library(purrr);
library(survival); library(doParallel)


#####select only the
source('../functionsForInteractionAnalysis.R')


cohort <- "CGGA693"; cohort.names <- "cgga693"; sampleType = "mut"
method <- "JADE"; dims = "10"


cross = T

if (cross == T) {
	load(glue("ICA-JADE-{dims}.crossCellTypes.{tolower(cohort)}.IDH{sampleType}-filtered-noNorm-optimalRank.Rda"))
}
if (cross == F) {
	load(glue("ICA-{method}-{dims}.cellTypes.{tolower(cohort)}.IDH{sampleType}-filtered-noNorm-optimalRank.Rda"))
}



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
cellTypes <- names(ica.res.list)

####parallel solution adapted from https://stackoverflow.com/questions/38318139/run-a-for-loop-in-parallel-in-r#######
#symmetricBins <- c(4,7,8)

cellTypeInteractionList <- list()
covariateType <- "technical" ## or "biological"
#covariateType <- "biological" ## or "technical"

for (cellType in cellTypes) {
	
	componentMatix1 <- ica.res.list[[cellType]]$A ###mixing matrix
	ncomp1 <- ncol(componentMatix1)
	colnames(componentMatix1) <- paste("Dim.", 1:ncomp1, sep = "")
	
	componentInteractionList <- NA
	
	if (cohort == "TCGA") {
		allSampleNames <- gsub("\\.", "-", rownames(componentMatix1))
		rownames(componentMatix1) <- allSampleNames
		desiredSamples <- TCGAquery_SampleTypes(allSampleNames, typesample = c("TP", "TR", "TM"))
	} else {
		if (grep("CGGA", cohort)) {
			allSampleNames <- rownames(componentMatix1)
			desiredSamples <- grep("CGGA", allSampleNames, value = T)
		}
	}

	componentMatix1 <- componentMatix1[rownames(componentMatix1) %in% desiredSamples, 1:ncomp1]
		
	if (cohort == "TCGA") { 
		rownames(componentMatix1) <- rownames(componentMatix1) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
	}
		
	lv.1 <- colnames(componentMatix1);
	
	sizeCount <- length(lv.1)
	cores=detectCores()
	cl <- makeCluster(cores[1] - 6) #not to overload your computer
	#cl <- makeCluster(9)
	registerDoParallel(cl)
	#for (componentIndex in 1:ncol(componentCombinations)) {
	componentInteractionList <- NULL
	componentInteractionList <- foreach(componentIndex=1:sizeCount, .packages=c('survival', 'dplyr', 'magrittr')) %dopar% 
			{
				comp1 <- componentMatix1[ ,componentIndex] #%>% set_names(samplenames)
				componentvalues <- data.frame(samples = names(comp1), comp1 = comp1)
			
				clinicalData <- merge(clinical, componentvalues, by = 'samples')
				tempList <- list()
				hr = NA; pvalue = NA
				if (covariateType == "biological") {
					if (sampleType == "IDHmut") {
						tryCatch({
							cox.out = coxph(Surv(time,status) ~ comp1 + age + sex + codel + grade, data = clinicalData)
							coefficients = summary(cox.out)$coefficients
							hr <- coefficients['comp1', 1]
							pvalue <- coefficients['comp1', 5]
					    },error=function(x){})
							
					} else {
						tryCatch({
							cox.out = coxph(Surv(time,status) ~ comp1 + age + sex + grade, data = clinicalData)
							coefficients = summary(cox.out)$coefficients
							hr <- coefficients['comp1', 1]
							pvalue <- coefficients['comp1', 5]
					    },error=function(x){})
					}
				}
				if (covariateType == "technical") {
					#message(type)
					tryCatch({
						cox.out = coxph(Surv(time,status) ~ comp1 + age + sex, data = clinicalData)
						coefficients = summary(cox.out)$coefficients
						hr <- coefficients['comp1', 1]
						pvalue <- coefficients['comp1', 5]
				    },error=function(x){})
				}
				tempList <- data.frame(HR = hr, p.value = pvalue)
				#nameInfo[[binIndex]] <- c(component1, component2, binIndex)
				tempList #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
					#}
			}
	#}
		#stop cluster
		stopCluster(cl)
		names(componentInteractionList) <- lv.1
		cellTypeInteractionList[[cellType]] <- componentInteractionList
		print(glue("done for {cellType}"))
}

component.interactions <- lapply(cellTypeInteractionList, function(x) bind_rows(x, .id = "Component"))
component.interactions <- lapply(component.interactions, function(x) {
	FDR = p.adjust(x$p.value, method = "fdr")
	x$FDR = FDR
	x
})
component.interactions <- bind_rows(component.interactions, .id = "cellType")

if (cross == T) {
	write.table(component.interactions, glue("survivalAnalysis.ICs-{method}-{dims}-{sampleType}-cross-{covariateType}-cross.txt"), sep = "\t", quote = F, row.names = F)
}
if (cross == F) {
	write.table(component.interactions, glue("survivalAnalysis.ICs-{method}-{dims}-{sampleType}-{cohort}-{covariateType}-{cohort}.txt"), sep = "\t", quote = F, row.names = F)
}

