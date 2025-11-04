library(argparse)

# Set up argument parser
parser <- ArgumentParser(description = "Perform regression analysis")
parser$add_argument("--ica_model_directory", type = "character", help = "a directory with .Rda files of the fitted ICA models", required = TRUE)
parser$add_argument("--is_tcga", action = "store_true", help = "use if the TCGA datasets are being analyzed", default = FALSE)
parser$add_argument("--clinical_file", type = "character", help = "a file containing the clinical data", required = TRUE)
parser$add_argument("--interaction_file", type = "character", help = "a file with all pairwise comparisons created by CSI-TME", required = TRUE)
parser$add_argument("--covariates", type = "character", help = "a comma separated list of covariates (default age,sex)", required = TRUE)
parser$add_argument("--strata", type = "character", help = "at most one variable for stratifing the model (optional)", required = F, default = NULL)
parser$add_argument("--output_file", type = "character", help = "a file to store the output file", required = T) #####default doesnt store the fitted models


args <- parser$parse_args()
ica_model_directory <- args$ica_model_directory
is_tcga <- args$is_tcga
covariates <- args$covariates
strata <- args$strata
clinical_file <- args$clinical_file
interaction_file <- args$interaction_file
output_file <- args$output_file


#library(clusterProfiler)

suppressPackageStartupMessages(
	{
		library(dplyr); library(magrittr);
		library(reshape2); library(survival); 
		library(glue); library(TCGAutils);
	}
)



source('codebase/R/functionsForInteractionAnalysis.R')


progress_monitor <- function(n_iter, current_iter, percent_offset = 5) {
	percent_complete <- (current_iter / n_iter) * 100
  
  	last_printed <- floor((current_iter - 1) / (n_iter * percent_offset / 100)) * percent_offset
	next_print <- floor(percent_complete / percent_offset) * percent_offset
  
	if (next_print > last_printed) {
		cat(sprintf("Progress: %d%% completed\n", next_print))
		flush.console()
	}
}


#####select only the

N_BOOTSRAPS = 10 ### curretnly we fix with 10 bootstraps. No need to change this unless you have specific reason

clinical_data <- read.table(clinical_file, sep = "\t", header = T)
ica_models_files <- list.files(ica_model_directory, full.name = T)

covariates <- strsplit(covariates, split = ",") %>% unlist()
covariates_str <- paste(covariates, collapse = " + ")

if (is.null(strata)) {
	cox_string <- paste("Surv(time, status)", " ~ ", "bins + comp1 + comp2 + ", covariates_str, sep = "")
} else {
	#cox_string <- paste("Surv(time, status)", " ~ ", "bins + comp1 + comp2 + ", covariates_str, "+ strata(strata)", sep = "")
	cox_string <- paste("Surv(time, status)", " ~ ", "bins + comp1 + comp2 + ", covariates_str, "+", " strata(", strata, ")", sep = "")
}
#cox_string <- paste("Surv(time, status)", " ~ ", "bins + comp1 + comp2 + ", covariates_str, sep = "")
cox_formula <- as.formula(cox_string)
#message(cox_formula)
		
ica.res.list <- list()
for (model_file in ica_models_files) {
	cell_type <- gsub("\\.Rda", "", basename(model_file))
	load(model_file)
	ica.res.list[[cell_type]] <- ica_model
}

binMapList <- componentBinningFunction(pcaList = ica.res.list, factorization = "ICA", sample.names = NULL)
cellTypeCombinations <- combn(names(binMapList), 2)


interaction_data <- read.table(interaction_file, sep = "\t", header = T)
interaction_data <- na.omit(interaction_data)
interaction_data <- interaction_data[interaction_data$FDR < 0.20, ]
interaction_data$type <- ifelse(interaction_data$HR > 0, "Pro-Tumor", "Anti-Tumor")
#nteractionFile <- interaction_data[interaction_data$pvalue < 0.007, ]

iteration.stats.list <- list()
for (row.index in 1:nrow(interaction_data)) {
	int.line <- interaction_data[row.index, ]
	sign.int <- sign(int.line$HR)
	#if (int.line$FDR < 0.30) {
		cellTypes <- strsplit(int.line$cellTypes, split = ":") %>% unlist()
		components <- strsplit(int.line$components, split = ":") %>% unlist()
		cellType1 <- cellTypes[1]; component1 = components[1]
		cellType2 <- cellTypes[2]; component2 = components[2]
		int.bin <- int.line[ ,'bin']
		
		score.1 <- ica.res.list[[cellType1]]$A %>% set_colnames(paste("IC", 1:ncol(.), sep = "."))
		score.2 <- ica.res.list[[cellType2]]$A %>% set_colnames(paste("IC", 1:ncol(.), sep = "."))
		
		comp1.values <- score.1[ ,component1]; comp1.bins <- binMapList[[cellType1]][ ,component1]
		comp2.values <- score.2[ ,component2]; comp2.bins <- binMapList[[cellType2]][ ,component2]
		
		if (all(rownames(comp1.values) == rownames(comp2.values))) {
			componentvalues <- data.frame(samples = names(comp1.values), comp1 = comp1.values, comp2 = comp2.values) %>% set_rownames(NULL)
			binValues <- data.frame(samples = names(comp1.bins), bin1 = comp1.bins, bin2 = comp2.bins) %>% set_rownames(NULL)
			
			if (is_tcga) {
				componentvalues$samples = gsub("\\.", "-", componentvalues$samples) 
				binValues$samples = gsub("\\.", "-", binValues$samples)
			
				#desiredSamples <- TCGAquery_SampleTypes(allSampleNames, typesample = c("TP", "TR", "TM"))
				primary <- TCGAsampleSelect(componentvalues$samples, "01")
				recurrant <- TCGAsampleSelect(componentvalues$samples, "02")
				metastatic <- TCGAsampleSelect(componentvalues$samples, "06")
				
				desiredSamples <- componentvalues$samples[primary | recurrant | metastatic]
			} else {
				allSampleNames <- rownames(binMap1)
				#desiredSamples <- grep("CGGA", allSampleNames, value = T)
				desiredSamples <- allSampleNames
			}		
			
			componentvalues <- componentvalues[componentvalues$samples %in% desiredSamples, ]
			binValues <- binValues[binValues$samples %in% desiredSamples, ]
						
			
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
			variables$samples <- TCGAutils::TCGAbarcode(variables$samples, participant = TRUE)
			indexes <- which(!duplicated(variables$samples))
			variables <- variables[indexes, ]
			
			clinicalData <- merge(clinical_data, variables, by = "samples")
			clinicalData <- clinicalData[!duplicated(clinicalData$samples), ]
			clinicalData[[strata]] <- factor(clinicalData[[strata]])
			iter.values <- list();
			for (iter in 1:(N_BOOTSRAPS + 1)) {
				if (iter == 1) {
					indexes <- sample(1:nrow(clinicalData), nrow(clinicalData), replace = F)
					dt1 <- clinicalData[indexes, ]
				} else {
					size <- floor(nrow(clinicalData) * 80 / 100)
					indexes <- sample(1:nrow(clinicalData), size, replace = F)
					dt1 <- clinicalData[indexes, ]
				}
				tryCatch({
					cox.out = coxph(cox_formula, data = dt1)
					coefficients = summary(cox.out)$coefficients
					hr <- coefficients['bins', 1]
					pvalue <- coefficients['bins', 5]
					values <- data.frame(hr = hr, pvalue = pvalue)
					iter.values[[iter]] <- values 
			    },error=function(x){})
			}
			iter.values <- bind_rows(iter.values, .id = "iteration")
			iter.stats <- sum(sign(iter.values$hr) == sign.int & iter.values$pvalue < 0.05)
			crossvalidation_accuracy <- round((iter.stats / N_BOOTSRAPS * 100), 2)
			#iteration.stats.list[[row.index]] <- iter.stats
			iteration.stats.list[[row.index]] <- crossvalidation_accuracy
		}
		#}
	progress_monitor(n_iter = nrow(interaction_data), current_iter = row.index, percent_offset = 5)
	#message(glue("done for {row.index} interactions"))
}

interaction_data$crossvalidation_accuracy <- unlist(iteration.stats.list)
#filtered.interactions <- interaction_data[interaction_data$Reproduciblity >= 8, ]
write.table(interaction_data, file = output_file, sep = "\t", quote = F, row.names = F)


