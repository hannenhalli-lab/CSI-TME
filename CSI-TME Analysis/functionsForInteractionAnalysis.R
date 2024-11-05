
get.bin.intervals <- function(cp,nbins,recursion = 1) { ##From Assafg
  binp = seq(from = 0, to = 1, by = 1/nbins)
  #binp <- c(0.00, 0.20, 0.60, 1.00)
  b   = unique(quantile(cp,binp,na.rm=TRUE))
  if(length(b[which(!is.na(b))])==0) return(NA)
  add = 0
  if( b[1] == 0 ) add = 1
  b = b[b!=0]
  if( add == 1 ) b = c(0,b)
  if( length(b) < nbins-1 ) return(NA)
  bin = b[1:(length(b)-1)]
  
  intervals = findInterval(cp, bin)-1
  counts    = table(intervals)

  len = length(cp)
  if( max(abs(counts-(len/nbins))) > max(1,0.005*len)
          && recursion > 0 ) {
    noise = sd(cp)/200
    cp = cp + runif(len, -noise, noise)
    return(get.bin.intervals(cp,nbins,recursion-1))
  }
  return(intervals)
}


componentBinningFunction <- function(pcaList, factorization, sample.names = NULL) {
	outputList <- list()
	if (factorization == "PCA") {
		for (index in 1:length(pcaList)) {
			fileName = pcaList %>% names() %>% .[index]
			cellType = gsub(".+Job\\d+_|_Window.+", "", fileName)
			pcaRes <- pcaList[[index]]
			componentMatix <- pcaRes$ind$coord
			if (ncol(componentMatix) >= 100) {
				#componentMatix <- componentMatix[ ,1:100]
				mappedBins <- apply(componentMatix, 2, get.bin.intervals, nbins = 3)
				outputList[[fileName]] <- mappedBins %>% set_rownames(rownames(componentMatix))
			} else {
				outputList[[fileName]] <- NA
			}
		}
	}
	if (factorization == "NMF") {
		for (index in 1:length(pcaList)) {
			fileName = pcaList %>% names() %>% .[index]
			cellType = gsub(".+Job\\d+_|_Window.+", "", fileName)
			pcaRes <- pcaList[[index]]
			componentMatix <- pcaRes@fit@H %>% t()
			if (ncol(componentMatix) >= 1) {
				#componentMatix <- componentMatix[ ,1:100]
				mappedBins <- apply(componentMatix, 2, get.bin.intervals, nbins = 3)
				col.names <- paste("Dim.", 1:ncol(mappedBins), sep = "")
				if (is.null(sample.names)) {
				  outputList[[fileName]] <- mappedBins %>% set_rownames(rownames(componentMatix)) %>% set_colnames(col.names)
				} else {
				  outputList[[fileName]] <- mappedBins %>% set_rownames(sample.names) %>% set_colnames(col.names)
				}
				
			} else {
				outputList[[fileName]] <- NA
			}
		}
	}
	if (factorization == "ICA") {
	  for (index in 1:length(pcaList)) {
	    fileName = pcaList %>% names() %>% .[index]
	    cellType = gsub(".+Job\\d+_|_Window.+", "", fileName)
	    pcaRes <- pcaList[[index]]
	    componentMatix <- pcaRes$A
	    if (ncol(componentMatix) >= 1) {
	      #componentMatix <- componentMatix[ ,1:100]
	      mappedBins <- apply(componentMatix, 2, get.bin.intervals, nbins = 3)
	      col.names <- paste("Dim.", 1:ncol(mappedBins), sep = "")
	      if (is.null(sample.names)) {
	        outputList[[fileName]] <- mappedBins %>% set_rownames(rownames(componentMatix)) %>% set_colnames(col.names)
	      } else {
	        outputList[[fileName]] <- mappedBins %>% set_rownames(sample.names) %>% set_colnames(col.names)
	      }
	      
	    } else {
	      outputList[[fileName]] <- NA
	    }
	  }
	}
	return(outputList)
}

projectionBinningFunction <- function(projectionList) {
	outputList <- list()
	for (index in 1:length(projectionList)) {
		fileName = projectionList %>% names() %>% .[index]
		cellType = gsub(".+Job\\d+_|_Window.+", "", fileName)
		componentMatix <- projectionList[[index]]
		if (!is.null(ncol(componentMatix))){
			if (ncol(componentMatix) >= 5) {
				#componentMatix <- componentMatix[ ,1:100]
				mappedBins <- apply(componentMatix, 2, get.bin.intervals, nbins = 3)
				rownames(mappedBins) <- rownames(componentMatix)
				outputList[[fileName]] <- mappedBins
			} else {
				outputList[[fileName]] <- NA
			}
		}  
	}
	return(outputList)
}

medianExpressionBinningFunction <- function(usedGEP, topGeneList) {
	binnedExpressionList <- list()
	medianExpressionList <- list()
	commonGenesList <- list()
	for (index in 1:length(usedGEP)) {
		fileName = usedGEP %>% names() %>% .[index]
		cellType = gsub(".+Job\\d+_|_Window.+", "", fileName)
		GEP <- usedGEP[[fileName]]
		topGenes <- topGeneList[[fileName]]
		if (ncol(GEP) >= 100) {
			tempMedianExpression <- list()
			tempCommonGenes <- list()
			for (i in 1:length(topGenes)) {
				topGeneNames <- topGenes[[i]]
				commonGenes <- GEP %>% .[ ,colnames(.) %in% topGeneNames] %>% colnames
				if (is.null(commonGenes)) {
					tempCommonGenes[[i]] <- NA
					tempMedianExpression[[i]] <- NA
				} else {
					tempCommonGenes[[i]] <- commonGenes
					medianExpression <- GEP %>% .[ ,colnames(.) %in% topGeneNames] %>% apply(., 1, median, na.rm = T)
					tempMedianExpression[[i]] <- medianExpression
				}
			}
			medianExpressionDf <- do.call("cbind", tempMedianExpression) %>% set_colnames(paste("Dim", rep(1:50), sep = "."))
			mappedBins <- apply(medianExpressionDf, 2, get.bin.intervals, nbins = 3)
			if (is.list(mappedBins)) {
				mappedBins <- do.call("cbind", mappedBins) %>% set_colnames(paste("Dim", rep(1:50), sep = "."))
			}
			binnedExpressionList[[fileName]] <- mappedBins
			medianExpressionList[[fileName]] <- medianExpressionDf
			commonGenesList[[fileName]] <- tempCommonGenes
			
		} else {
			binnedExpressionList[[fileName]] <- NA
			medianExpressionList[[fileName]] <- NA
			commonGenesList[[fileName]] <- NA
		}
	}
	return(list(commonGenes = commonGenesList, medianExpression = medianExpressionList, binnedExpression = binnedExpressionList))
}

####used from Assaf's
getInteraction <- function (tt) {
	bins.binary = sapply(1:9, function(tt) ifelse(bins == tt, 1, 0))
	colnames(bins.binary) = sapply(1:9,function(v) paste('bin',v,sep=''))
}

interaction.membership.mapping <- function(binningList, cellType1, component1, cellType2, component2, intBin) {
	
	binMap1 <- binningList[[cellType1]] #%>% set_rownames(samplenames1)
	binMap2 <- binningList[[cellType2]] #%>% set_rownames(samplenames1)
	
	if (all(rownames(binMap1) == rownames(binMap2))) {
		
		if (cohort == "TCGA") {
			allSampleNames <- gsub("\\.", "-", rownames(binMap1))
			rownames(binMap1) <- allSampleNames
			rownames(binMap2) <- allSampleNames
			
			desiredSamples <- TCGAquery_SampleTypes(allSampleNames, typesample = c("TP", "TR", "TM"))
		} else {
			if (length(grep("CGGA", cohort)) == 1) {
				allSampleNames <- rownames(binMap1)
				desiredSamples <- grep("CGGA", allSampleNames, value = T)
			}
		}
	}

	binMap1 <- binMap1[rownames(binMap1) %in% desiredSamples, ]
	binMap2 <- binMap2[rownames(binMap2) %in% desiredSamples, ] 
	
	
	if (cohort == "TCGA") { 
		rownames(binMap1) <- rownames(binMap1) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
		rownames(binMap2) <- rownames(binMap2) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
	}
	
	
	bin.map <- cbind(binMap1[ ,component1], binMap2[ ,component2])
	
	if (intBin == "Bin.1") {
		prop.comp1 <- sum(bin.map[ ,1] == 0) / nrow(bin.map)
		prop.comp2 <- sum(bin.map[ ,2] == 0) / nrow(bin.map)
	}
	
	if (intBin == "Bin.3") {
		prop.comp1 <- sum(bin.map[ ,1] == 2) / nrow(bin.map)
		prop.comp2 <- sum(bin.map[ ,2] == 0) / nrow(bin.map)
	}
	
	if (intBin == "Bin.7") {
		prop.comp1 <- sum(bin.map[ ,1] == 0) / nrow(bin.map)
		prop.comp2 <- sum(bin.map[ ,2] == 2) / nrow(bin.map)
	}
	
	if (intBin == "Bin.9") {
		prop.comp1 <- sum(bin.map[ ,1] == 2) / nrow(bin.map)
		prop.comp2 <- sum(bin.map[ ,2] == 2) / nrow(bin.map)
	}
	
	expectation <- prop.comp1 * prop.comp1
	
    #bins = bin.map[ ,1] * 3 + bin.map[ ,2] + 1 #####this one is used in spageFinder
	bins = bin.map[ ,2] * 3 + bin.map[ ,1] + 1 ####I switched the indexes for calculation to call bins in order of comp1:comp2
  	bins.binary = sapply(1:9, function(tt) ifelse(bins == tt, 1, 0))
    colnames(bins.binary) = sapply(1:9,function(v) paste('Bin',v,sep='.'))
	bins.binary <- data.frame(bins.binary) %>% cbind(samples = rownames(bins.binary), .)
	clinicalData = merge(clinical, bins.binary, by.x = "samples")

	interactionStats <- clinicalData[ ,intBin] %>% set_names(clinicalData$samples)
	interactingSamples <- interactionStats[interactionStats == 1] %>% names
	
	return(list(membership = interactionStats, members = interactingSamples, expectation = expectation))
	
}

perform.cox <- function(binaryBins,index, clinicalData, components = T, type, covariate) {
	bins <- binaryBins[ ,index] %>% data.frame(samples = names(.), bins = .)
	#dt1 = cbind(clinical, bin = bins)
	dt1 <- merge(clinicalData, bins, by.x = "samples", by.y = "samples") %>% unique()
	indexes <- which(!duplicated(dt1$samples))
	dt1 <- dt1[indexes, ]
		#signed.delta.loglik = NA
		hr = NA; pvalue = NA
		if (components) {
			tryCatch({
				#cox.out = coxph(Surv(time,status) ~ bins + comp1 + comp1 + age + strata(type), data = dt1)
				if (covariate == "biological") {
					if (type == "IDHmut") {
						cox.out = coxph(Surv(time,status) ~ bins + comp1 + comp2 + age + sex + codel + grade, data = dt1)
					} else {
						cox.out = coxph(Surv(time,status) ~ bins + comp1 + comp2 + age + sex + grade, data = dt1)
					}
				}
				if (covariate == "technical") {
					#message(type)
					cox.out = coxph(Surv(time,status) ~ bins + comp1 + comp2 + age + sex, data = dt1)
				}
				coefficients = summary(cox.out)$coefficients
				hr <- coefficients['bins', 1]
				pvalue <- coefficients['bins', 5]
		    },error=function(x){})
		} else {
			tryCatch({
				cox.out = coxph(Surv(time,status) ~ bins + age, data = dt1)
				coefficients = summary(cox.out)$coefficients
				hr <- coefficients['bins', 1]
				pvalue <- coefficients['bins', 5]
		    },error=function(x){})
		}
	    options(warn=1)
	    return(c(hr, pvalue))
  }

sampleNameExtractor <- function(cellType, gepList) {
samplenames <- rownames(gepList[[cellType]]) %>% 
		gsub("\\.", "-", .) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
		
return(samplenames)
}


outputParser <- function(inputList, parameterList) {	
	#componentLength = parameterList[['componentLength']]
	selection = parameterList[['selection']] 
	componentNames = parameterList[['componentNames']]
	binNames = parameterList[['binNames']]
	
	if (all(is.na(inputList))) {
		
	}else {
		componentLength <- length(inputList[[1]])
		componentNames <- names(inputList[[1]])
		fun <- function (index, data, selection) {
			dataList <- data[[1]][[index]] %>% unlist() %>% matrix(., ncol = 2, byrow = T) #%>% data.frame()
			dataList[,selection]
		}
		#tmpList <- inputList
		sapply(1:componentLength, fun, selection = selection, data = inputList) %>% t() %>%
					data.frame() %>% set_rownames(componentNames) %>% set_colnames(binNames)
	}	
}



symmetricBinCorrection <- function (InputData, select, mapBins) {
	symmetricBins = c("Bin.4", "Bin.7", "Bin.8")
	conversion <- c(Bin.4 = "Bin.2", Bin.7 = "Bin.3", Bin.8 = "Bin.6")
	indexes <- InputData$bin %in% symmetricBins
	if (sum(indexes) > 0) {
		if (select == "components") {
			componentPairs <- InputData[indexes, "components"] %>% strsplit(., ":") %>% 
						unlist() %>% matrix(., ncol = 2, byrow = T) %>% data.frame()
			componentPairs <- paste(componentPairs[,2], componentPairs[,1], sep = ":")
			InputData[indexes, 'components'] <- componentPairs	
		}
		if (select == "cellTypes") {
			cellTypePairs <- InputData[indexes, "cellTypes"] %>% strsplit(., ":") %>% 
						unlist() %>% matrix(., ncol = 2, byrow = T) %>% data.frame()
			cellTypePairs <- paste(cellTypePairs[,2], cellTypePairs[,1], sep = ":")
			InputData[indexes, 'cellTypes'] <- cellTypePairs
		
		}
		if (mapBins) {
			binNames <- InputData[indexes, 'bin'] %>% as.character()
			InputData[indexes, 'bin'] <- conversion[binNames]
		}
		return(InputData)
	} else {
		return(InputData)
	}
	
}

fdrFunction <- function(x, method) {
	tempData <- x[ ,c(1,2)]
	bin.names <- colnames(x)[3:ncol(x)]
	fdrData <- x[ ,-c(1,2)] %>% data.frame()#%>% dim
	fdr <- apply(fdrData, 2, p.adjust, method) %>% data.frame() %>% set_colnames(bin.names) %>%
		 mutate(cellTypes = tempData[ ,1], components = tempData[ ,2], .before = bin.names[1])
	
	rownames(fdr) <- paste(fdr[,1], fdr[ ,2],  sep = "") ####to avoid duplicate rowname error while unsplitting the data
	return(fdr)
}

getSignificantInteractions <- function(hazardRatiosDf, pvalueDf, method, threshold, factorList, keep.bins = c("Bin.1", "Bin.3", "Bin.7", "Bin.9"), symmetricBinCorrection = T) {
	#factorList <- hazardRatiosDf[ ,factorList] %>% list(.[,1], .[,2])
	
	hazardRatiosDf <- hazardRatiosDf[ ,colnames(hazardRatiosDf) %in% c('cellTypes', 'components' ,keep.bins)]
	pvalueDf <- pvalueDf[ ,colnames(pvalueDf) %in% c('cellTypes', 'components' ,keep.bins)]
	
	meltedHazardsDf <- reshape2::melt(hazardRatiosDf, id.vars = c("cellTypes", "components")) %>% set_colnames(c("cellTypes", "components", "bin", "HR"))
	
	if (symmetricBinCorrection == T) {
		meltedHazardsDf <- symmetricBinCorrection(InputData = meltedHazardsDf, select = "components", mapBins = F) 
		meltedHazardsDf <- symmetricBinCorrection(InputData = meltedHazardsDf, select = "cellTypes",  mapBins = T) ###mapBins argument is used only once and replace all the symmetric bin names
		#meltedHazardsDf %>% split(. , .[ ,factorList]) %>% lapply
	}
	
	#fdrDf <- pvalueDf[ ,-c(1,2)] %>% apply(., 2, p.adjust, method) %>% data.frame() %>%
	#		mutate(cellTypes = pvalueDf[ ,1], components = pvalueDf[ ,2], .before = "Bin.1")
	factorList <- pvalueDf[[factorList]] %>% as.factor
	fdrDf <- pvalueDf %>% split(., factorList) %>% lapply(., function(x) fdrFunction(x, method = method)) %>% unsplit(., factorList)
	rownames(fdrDf) <- NULL 		
		
	meltedFdrDf <- reshape2::melt(fdrDf, id.vars = c("cellTypes", "components")) %>% set_colnames(c("cellTypes", "components", "bin", "FDR"))
	
	if (symmetricBinCorrection == T) {
		meltedFdrDf <- symmetricBinCorrection(InputData = meltedFdrDf, select = "components", mapBins = F)
	}
	
	meltedHazardsDf$FDR <- meltedFdrDf$FDR
	meltedHazardsDf$pvalue <- reshape2::melt(pvalueDf)$value
	#significantInterations <- fdrDf[ ,-c(1,2)] %>% apply(., 2, function(x) x < threshold) %>% melt()
	#meltedHazardsDf <- meltedHazardsDf[significantInterations$value, ]
	#meltedFdrDf <- meltedFdrDf[significantInterations$value, ]
	#meltedHazardsDf$FDR <- meltedFdrDf$value
	significantInterations <- meltedHazardsDf[meltedHazardsDf$FDR < threshold, ]
	return(significantInterations)
}

getTopGenes <- function(decomposition.list, analysisType = "PLIER") {
	if (analysisType == "PLIER") {
		lodings <- decomposition.list$Z
		componenetNames <- paste("LV", 1:ncol(lodings), sep = "")
		sortedSites <- apply(lodings, 2, function(x) {
			y = x[order(x, decreasing = T)]; 
			nonzero = sum(x != 0); 
			y[1:nonzero] %>% names
		}) 
	}
	if (analysisType == "PCA") {
		lodings <- sweep(decomposition.list$var$coord,2,sqrt(decomposition.list$eig[1:ncol(decomposition.list$var$coord),1]),FUN="/")
		componenetNames <- paste("Dim", 1:ncol(lodings), sep = ".")
		sortedSites <- apply(lodings, 2, function(x) {
			y = x[order(x, decreasing = T)]; 
			nonzero = sum(x != 0); 
			y[1:nonzero] %>% names
		}) 
	}
	sortedSites <- sortedSites %>% set_colnames(componenetNames)
	return(sortedSites)
}

mutualInformation <- function(binningList, cellType1, component1, cellType2, component2, intBin) {
	binMap1 <- binningList[[cellType1]] #%>% set_rownames(samplenames1)
	binMap2 <- binningList[[cellType2]] #%>% set_rownames(samplenames1)
	
	if (all(rownames(binMap1) == rownames(binMap2))) {
		
		if (cohort == "TCGA") {
			allSampleNames <- gsub("\\.", "-", rownames(binMap1))
			rownames(binMap1) <- allSampleNames
			rownames(binMap2) <- allSampleNames
			
			desiredSamples <- TCGAquery_SampleTypes(allSampleNames, typesample = c("TP", "TR", "TM"))
		} else {
			if (grep("CGGA", cohort)) {
				allSampleNames <- rownames(binMap1)
				desiredSamples <- grep("CGGA", allSampleNames, value = T)
			}
		}
		
	}

	binMap1 <- binMap1[rownames(binMap1) %in% desiredSamples, ]
	binMap2 <- binMap2[rownames(binMap2) %in% desiredSamples, ] 
	
	
	if (cohort == "TCGA") { 
		rownames(binMap1) <- rownames(binMap1) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
		rownames(binMap2) <- rownames(binMap2) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
	}
	
	
	component1 <- binMap1[ ,component1]
	component2 <- binMap2[ ,component2]
	
	
	if (intBin == "Bin.1") {
		vec1 <- ifelse(component1 == 0, 1, 0)
		vec2 <- ifelse(component2 == 0, 1, 0)
	}
	if (intBin == "Bin.3") {
		vec1 <- ifelse(component1 == 2, 1, 0)
		vec2 <- ifelse(component2 == 0, 1, 0)
	}
	if (intBin == "Bin.7") {
		vec1 <- ifelse(component1 == 0, 1, 0)
		vec2 <- ifelse(component2 == 2, 1, 0)
	}
	if (intBin == "Bin.9") {
		vec1 <- ifelse(component1 == 2, 1, 0)
		vec2 <- ifelse(component2 == 2, 1, 0)
	}
	mut.information <- infotheo::mutinformation(vec1, vec2, method = "emp")
	return(mut.information)
}

get.LR.interactions <- function(list1, list2, database) {
	lr.total <- nrow(database)
	database <- split(database, database$Ligand)
	
	ligands.all <- names(database) %>% unique
	receptor.all <- lapply(database, function(x) x$Receptor) %>% unlist() %>% unique()
	
	ligands.total <- length(ligands.all)
	receptor.total <- length(receptor.all)
	
	length.comp1 <- length(list1)
	length.comp2 <- length(list2)
	
	
	
	ligands.comp1 <- list1[list1 %in% ligands.all] #####ligands in component 1
	receptors.comp2 <- list2[list2 %in% receptor.all]  #####receptors in component 2
	receptors <- database[ligands.comp1] %>% lapply(., function(x) x$Receptor) #####receptor of ligands in component 1
	
	
	ligands.comp1.length <- length(ligands.comp1) ####total number of unique ligands in component 1
	receptors.comp2.length <- length(receptors.comp2) ####total number of unique receptors in component 2
	
	
	LR <- lapply(receptors, function(x) x[x%in% list2] %>% data.frame(Gene2 = .)) #####get receptors in component 2
	LR <- bind_rows(LR, .id = "Gene1") %>% unique()####make LR pairs
	
	
	if (nrow(LR) == 0) {
		LR <- data.frame(Gene1 = NA, Gene2 = NA)
	}
	LR$type <- "LR" 
	LR$size.comp1 <- length.comp1; LR$size.comp2 <- length.comp2
	LR$size.ligands <- ligands.comp1.length; LR$size.receptors <- receptors.comp2.length
	
	ligands.comp2 <- list2[list2 %in% ligands.all] #####ligands of component 2
	receptors.comp1 <- list1[list1 %in% receptor.all]  #####receptors in component 1
	receptors <- database[ligands.comp2] %>% lapply(., function(x) x$Receptor) #####receptor of ligands in component 2
	
	ligands.comp2.length <- length(ligands.comp2) ####total number of unique ligands in component 1
	receptors.comp1.length <- length(receptors.comp1) ####total number of unique receptors in component 2
	
	
	RL <- lapply(receptors, function(x) x[x%in% list1] %>% data.frame(Gene2 = .)) #####get receptors in component 1
	RL <- bind_rows(RL, .id = "Gene1") %>% unique() ####make LR pairs
	if (nrow(RL) == 0) {
		RL <- data.frame(Gene1 = NA, Gene2 = NA)
	}
	RL <- RL[ ,c(2,1)] %>% set_names(c("Gene1", "Gene2")) ### switch indexes to make RL pairs
	RL$type <- "RL"
	RL$size.comp1 <- length.comp1; RL$size.comp2 <- length.comp2
	RL$size.ligands <- ligands.comp2.length; RL$size.receptors <- receptors.comp1.length
	
	rbind(LR, RL)
}

format.interactions <- function(interactions.df, input = "cellTypePairs") {
	if (input == "cellTypePairs") {
		tidyr::separate(interactions.df, col = cellTypes, into = c("cellType1", "cellType2"), sep = ":|-") %>% 
			tidyr::separate(., col = components, into = c("component1", "component2"), sep = ":|-") %>%
			mutate(., partner1 = paste(cellType1, component1, sep = ":"), partner2 = paste(cellType2, component2, sep = ":"))
	}
}

get.partners <- function(inputDf) {
	###must be outputed from format.interactions format####
	partners.df <- inputDf[ ,c('partner1', 'partner2')] %>% apply(., 1, sort) %>% t() %>% data.frame()
	colnames(partners.df) <- c('partner1', 'partner2')
	partners.df$partners <- paste(partners.df$partner1, partners.df$partner2, sep = "_")
	inputDf$partners <- partners.df$partners
	return(inputDf)
}

####do a fisher test####
fisherTestFunction <- function(list1, list2, universe) {
	
	if (is.numeric(universe)) {
		universe = universe
	} else {
		list1 <- list1[list1 %in% universe]
		list2 <- list2[list2 %in% universe]
		universe <- length(universe)
	}
	
	c1 <- sum(list1 %in% list2)
	c2 <- length(list1) - c1
	c3 <- length(list2) - c1
	c4 <- universe - (c3 + length(list1))
	
	fisherMatrix <- matrix(c(c1,c2,c3,c4), nrow = 2)
	fisherTest <- fisher.test(fisherMatrix)
	enrichment <- fisherTest$estimate
	pvalue <- fisherTest$p.value
	CI <- fisherTest$conf.int
	
	return(c(L1 = length(list1), L2 = length(list2), common = c1, enrichment = enrichment, pvalue = pvalue, CI = CI))
}