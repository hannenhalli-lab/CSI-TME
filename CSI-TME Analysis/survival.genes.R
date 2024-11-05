#library(clusterProfiler)
library(dplyr); library(magrittr);
library(gtools); library(preprocessCore); 
library(reshape2); library(ggplot2);
library(glue); library(TCGAbiolinks); library(purrr);
library(survival); library(doParallel)
library(preprocessCore)


#####select only the
source('../functionsForInteractionAnalysis.R')

sampleType = "mut"
cohort <- "CGGA693"; 
expression.path <- glue("{tolower(cohort)}_{sampleType}_filtered/")
expression.files <- c("Bcell.txt","Endothelial.txt","Malignant.txt","Myeloid.txt","Oligos.txt","Stromal.txt","Tcell.txt")

if (cohort == "TCGA") {
	cancerType <- c("GBM", "LGG")
	clinical <- read.table("/Users/singha30/DataAndAnalysis/SPAGE/RB1-Project/AveragedSamples/AllAvialableClinicalDataMod.txt", sep = "\t", header = T)
	clinical <- clinical[clinical$type %in% cancerType, ]
	
	codel <- read.csv("../All_metadata-CCGA.csv", header = T) %>% set_names(c('samples' ,'project', 'time', 'status', 'IDH', 'grade', 'age', 'subtype', 'codel'))
	codel <- codel[grep("TCGA", codel$samples), c('IDH', 'samples', 'codel', 'grade')]
	codel$samples <- gsub("\\.", "-", codel$samples) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
	clinical <- merge(clinical, codel, by = "samples")
}

if (cohort == "CGGA693") {
	#componentMatrix1 <- read.table(glue("../raw_gene_exp_matrix/cgga_693_idh_{sampleType}_rsem.txt"), header = T, row.names = 1)
	clinical <- read.csv("../All_metadata-CCGA.csv", header = T) %>% set_names(c('samples' ,'project', 'time', 'status', 'IDH', 'grade', 'age', 'subtype', 'codel'))
	#gender <- read.table("cgga.combinedMetadata.txt", sep = "\t", header = T) %>% .[ ,c('samples', 'Gender')]
	gender <- read.table("cgga.clinical.gender.txt", sep = "\t", header = T) %>% .[ ,c('samples', 'Gender')]
	clinical <- merge(clinical, gender, by = "samples")
	names(clinical)[names(clinical) == 'Gender'] <- "sex"
}

if (cohort == "CGGA325") {
	
	clinical <- read.csv("../All_metadata-CCGA.csv", header = T) %>% set_names(c('samples' ,'project', 'time', 'status', 'IDH', 'grade', 'age', 'subtype', 'codel'))
	gender <- read.table("cgga.combinedMetadata.txt", sep = "\t", header = T) %>% .[ ,c('samples', 'Gender')]
	clinical <- merge(clinical, gender, by = "samples")
	names(clinical)[names(clinical) == 'Gender'] <- "sex"
}




covariateType <- "technical" ## or "biological"
#covariateType <- "biological" ## or "technical"


survival.list <- list()
for (expression.file in expression.files) {
	#componentMatrix1 <- read.table(glue("../raw_gene_exp_matrix/tcga_idh_{sampleType}_tpm.txt"), header = T, row.names = 1)
	componentMatrix1 <- read.table(glue("{expression.path}/{expression.file}"), header = T, row.names = 1)
	componentMatrix1 <- log2(componentMatrix1 + 0.001)
	componentMatrix1 <- normalize.quantiles(data.matrix(componentMatrix1)) %>% 
		set_rownames(rownames(componentMatrix1)) %>% set_colnames(colnames(componentMatrix1))
	
	
	componentMatrix1 <- t(componentMatrix1) %>% data.frame()
	componentInteractionList <- NA
	
	if (cohort == "TCGA") {
		allSampleNames <- gsub("\\.", "-", rownames(componentMatrix1))
		rownames(componentMatrix1) <- allSampleNames
		desiredSamples <- TCGAquery_SampleTypes(allSampleNames, typesample = c("TP", "TR", "TM"))
	}
	if (cohort == "CGGA693") {
		allSampleNames <- rownames(componentMatrix1)
		desiredSamples <- grep("CGGA", allSampleNames, value = T)
	}
	if (cohort == "CGGA325") {
		allSampleNames <- rownames(componentMatrix1)
		desiredSamples <- grep("CGGA", allSampleNames, value = T)
	}
	
	componentMatrix1 <- componentMatrix1[rownames(componentMatrix1) %in% desiredSamples, ]
		
	if (cohort == "TCGA") { 
		sample.names <- rownames(componentMatrix1) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
		indexes <- !duplicated(sample.names)
		componentMatrix1 <- componentMatrix1[indexes, ]
		rownames(componentMatrix1) <- sample.names[indexes]
		
	}
		
	lv.1 <- colnames(componentMatrix1);
	
	sizeCount <- length(lv.1)
	cores=detectCores()
	cl <- makeCluster(cores[1] - 6) #not to overload your computer
	#cl <- makeCluster(9)
	registerDoParallel(cl)
	#for (componentIndex in 1:ncol(componentCombinations)) {
	componentInteractionList <- NULL
	componentInteractionList <- foreach(componentIndex=1:sizeCount, .packages=c('survival', 'dplyr', 'magrittr')) %dopar% 
			{
				comp1 <- componentMatrix1[ ,componentIndex] #%>% set_names(samplenames)
				componentvalues <- data.frame(samples = rownames(componentMatrix1), comp1 = comp1)
			
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
		
		component.interactions <- bind_rows(componentInteractionList, .id = "Component")
		component.interactions$FDR <- p.adjust(component.interactions$p.value, method = "fdr")
		
		survival.list[[expression.file]] <- component.interactions
		
		print(glue("done for {expression.file}"))
} 


survivalDf <- bind_rows(survival.list, .id = "cellType")
write.table(survivalDf, glue("survivalAnalysis.Genes-IDH{sampleType}-{cohort}-cellTypes.txt"), sep = "\t", quote = F, row.names = F)

#write.table(component.interactions, glue("survivalAnalysis.Genes-IDH{sampleType}-{cohort}-{covariateType}-{cohort}.txt"), sep = "\t", quote = F, row.names = F)


#write.table(component.interactions, glue("survivalAnalysis.Genes-IDH{sampleType}-{cohort}-{covariateType}-{expression.file}"), sep = "\t", quote = F, row.names = F)


############################
##########following modified from https://r-coder.com/correlation-plot-r/
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
	data.cor <- cbind(x, y) %>% na.omit() %>% data.frame()
    Cor <- abs(cor(data.cor$x, data.cor$y)) # Remove abs function if desired
    txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
    if(missing(cex.cor)) {
        cex.cor <- 0.4 / strwidth(txt)
    }
    text(0.5, 0.5, txt,
         cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}



#files <- list.files(path = ".", pattern = "survivalAnalysis.Genes-IDHmut") 
#files <- list.files(path = ".", pattern = "Malignant.txt") %>% grep("IDHwt", ., value = T)
files <- list.files(path = ".", pattern = "cellTypes.txt") %>% grep("IDHmut", ., value = T)
file.names <- gsub("survivalAnalysis.Genes-|-cellTypes.txt", "", files)
names(files) <- file.names

survival.data <- lapply(files, read.table, sep = "\t", header = T)

survival.data <- survival.data %>% purrr::reduce(left_join, by = c("cellType", "Component"))
ggplot(data = survival.data, aes(x = HR.x, y =  HR.y)) + geom_point() + facet_wrap(~cellType) + ggpubr::stat_cor()

survival.data.HR <- survival.data[ ,c('HR.x', 'HR.y', 'HR')] %>% set_names(file.names) 



	  
	  
svglite::svglite(file = glue("plots/HazardRatio.Genes-Cohorts-{sampleType}.svg"), width = 5, height = 5)
pairs(survival.data.HR, upper.panel = panel.cor,  lower.panel = panel.smooth)
dev.off()



##########survival.datasurvival.datasurvival.datasurvival.datasurvival.data

files <- list.files(path = ".", pattern = "cellTypes.txt") %>% grep("IDHmut", ., value = T)
file.names <- gsub("survivalAnalysis.Genes-|-cellTypes.txt", "", files)
names(files) <- file.names

survival.data <- lapply(files, read.table, sep = "\t", header = T)

cgga <- survival.data[[1]] %>% subset(., select = c(cellType, Component, HR)) %>% split(., .$cellType) %>% lapply(., function(x) x[ ,-1])
tcga <- survival.data[[2]] %>% subset(., select = c(cellType, Component, HR)) %>% split(., .$cellType) %>% lapply(., function(x) x[ ,-1])

cellTypeCombinations <- combn(names(tcga), 2)

self.list <- list()
for (cellType in names(tcga)) {
	tcga.data <- tcga[[cellType]]
	cgga.data <- cgga[[cellType]]
	correlation.data <- merge(tcga.data, cgga.data, by = "Component")
	correlation <- cor(correlation.data[ ,'HR.x'], correlation.data[ ,'HR.y'])
	self.list[[cellType]] <- correlation
}
cross.correlation <- list()
for (index in 1:ncol(cellTypeCombinations)) {
	cellType1 <- cellTypeCombinations[1, index]
	cellType2 <- cellTypeCombinations[2, index]
	
	tcga.data <- tcga[[cellType1]]
	cgga.data <- cgga[[cellType2]]
	correlation.data <- merge(tcga.data, cgga.data, by = "Component")
	correlation <- cor(correlation.data[ ,'HR.x'], correlation.data[ ,'HR.y'])
	cross.correlation[index] <- correlation
}

data.plot1 <- data.frame(comparison = "same cell type", correlation = unlist(self.list))
data.plot2 <- data.frame(comparison = "across cell type", correlation = unlist(cross.correlation))

data.plot <- rbind(data.plot1, data.plot2)
data.plot$comparison <- factor(data.plot$comparison, levels = unique(data.plot$comparison))

svglite::svglite(file = glue("plots/HazardRatio.correlation.genes-acrossCohorts-{sampleType}.svg"), width = 3.5, height = 3.5)
ggplot(data = data.plot, aes(x = comparison, y = correlation)) + geom_boxplot() + stat_compare_means() + theme_bw()
dev.off()


data.plot <- data.frame(tcga = melt(mat2)$value, cgga = melt(mat1)$value) %>% data.matrix()
data.plot <- reshape2::melt(data.plot)



data.cor.list <- list()
for (cellType in cellTypes) {
	discovery <- survival.data[['TCGA']] %>% .[grepl(cellType, .$cellType), ]
	validation <- survival.data[['cross']] %>% .[grepl(cellType, .$cellType), ]
	
	data.cor <- merge(discovery, validation, by = "Component")
	data.cor <- data.cor %>% group_by(cellType.y) %>% summarize(cor(HR.x, HR.y, method = "spearman")) %>% data.frame() %>% set_names(c("cellType", "correlation"))
	data.cor.list[[cellType]] <- data.cor
}




files <- list.files(path = ".", pattern = "survivalAnalysis.PCs-IDHmut") 
file.names <- gsub(".+-|.txt", "", files)
names(files) <- file.names

survival.data <- lapply(files, read.table, sep = "\t", header = T)

survival.data <- survival.data %>% purrr::reduce(left_join, by = c("cellTypes", "Component"))

survival.data.HR <- survival.data[ ,c('HR.x', 'HR.y', 'HR')] %>% set_names(file.names) 



svglite::svglite(file = glue("plots/HazardRatio.correlation.genesVsICs-{sampleType}.svg"), width = 3.5, height = 3.5)
ggplot(data = data.plot, aes(x = comparison, y = correlation)) + geom_boxplot() + stat_compare_means() + theme_bw()
dev.off()
