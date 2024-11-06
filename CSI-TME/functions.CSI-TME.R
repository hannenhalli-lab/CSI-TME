
get.bin.intervals <- function(cp,nbins,recursion = 1) { 
  ##From Assafg
  #https://github.com/asmagen/SPAGEfinder
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



####used from Assaf's
getInteraction <- function (tt) {
	bins.binary = sapply(1:9, function(tt) ifelse(bins == tt, 1, 0))
	colnames(bins.binary) = sapply(1:9,function(v) paste('bin',v,sep=''))
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
