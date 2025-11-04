#library(clusterProfiler)

library(argparse)

# Set up argument parser
parser <- ArgumentParser(description = "get interaction penetrance and interaciton load")
parser$add_argument("--discovery_cohort", type = "character", help = "name of the discovery cohort", required = TRUE)
parser$add_argument("--target_cohort", type = "character", help = "name of the target cohort", required = TRUE)
parser$add_argument("--is_cohort", type = "character", help = "whether the cohort is 'discovery' or 'validation'", required = T)
parser$add_argument("--sample_file", type = "character", help = "a complete path to file with two columns containing the classification for samples (samples, class)", default = NULL)


library(dplyr); library(magrittr);
source("codebase/R/functionsForInteractionAnalysis.R")

args <- parser$parse_args()
cohort_name <- args$cohort_name
target_cohort <- args$target_cohort
sample_file <- args$sample_file
is_cohort <- args$is_cohort

if (is_cohort == "discovery") {
	model_directory <- paste("output/", project, "/ica_models/", sep = "")
}
if (is_cohort == "validation") {
	model_directory <- paste("output/validation/", target_cohort, "/ica_models/", sep = "")
}

interaction_file <- paste("output/", project, "/CSI-TME_significant_crossvalidation.txt", sep = "")
output_load <- paste("output/", project, "/interaction_load_", target_cohort,  ".txt", sep = "")
output_penetrance <- paste("output/", project, "/interaction_penetrance_", target_cohort,  ".txt", sep = "")

model_files <- list.files(model_directory, full.names = T)
cell_types <- gsub(".+/|\\..+", "", model_files)

ica.res.list <- list()

for (i in 1:length(model_files)) {
	model_file <- model_files[i]
	cell_type <- cell_types[i]
	load(model_file)
	ica.res.list[[cell_type]] <- ica_model
}


binMapList <- componentBinningFunction(pcaList = ica.res.list, factorization = "ICA")
interaction_data <- read.table(interaction_file, sep = "\t", header = T)
interactions.robust <- interaction_data %>% subset(., FDR < 0.20 & crossvalidation_accuracy > 70)
interactions.robust$string <- paste(interactions.robust$cellTypes, interactions.robust$components, interactions.robust$bin, sep = "_")

if (!is.null(sample_file)) {
	sample_data <- read.table(sample_file, sep = "\t", header = T)
	colnames(sample_data) <- c("samples", "class")
} else {
	sample_data <- NULL
}

progressIndex = 0
sample.matrix.list <- list()
for (index in 1:nrow(interactions.robust)) {
	interaction.string <- interactions.robust[index, ]
	
	cellTypes <- strsplit(interaction.string$cellTypes, ":") %>% unlist()
	components <- strsplit(interaction.string$components, ":") %>% unlist()
	int.bin <- interaction.string$bin
	
	cellType1 <- cellTypes[1] %>% gsub(".csv", "", .)
	cellType2 <- cellTypes[2] %>% gsub(".csv", "", .)
	
	dim1 <- components[1]
	dim2 <- components[2]
	
	comp1.bins <- binMapList[[cellType1]][ ,dim1]
	comp2.bins <- binMapList[[cellType2]][ ,dim2]
	
	binValues <- data.frame(barcodes = names(comp1.bins), bin1 = comp1.bins, bin2 = comp2.bins) %>% set_rownames(NULL)
	
	if (int.bin == "Bin.1") {
		bins <- ifelse(binValues$bin1 == 0 & binValues$bin2 == 0, 1, 0)
	}
	if (int.bin == "Bin.3") {
		bins <- ifelse(binValues$bin1 == 2 & binValues$bin2 == 0, 1, 0)
	}
	if (int.bin == "Bin.9") {
		bins <- ifelse(binValues$bin1 == 2 & binValues$bin2 == 2, 1, 0)
	}
	sample.matrix.list[[index]] <- bins %>% set_names(names(comp1.bins)) %>% data.frame()
}


sample.matrix.df <- do.call("cbind", sample.matrix.list) %>% set_colnames(interactions.robust$string) %>% t()

sample.matrix.df <- sample.matrix.df %>% data.matrix %>% reshape2::melt()
sample.matrix.df$Var2 <- gsub("\\.", "-", sample.matrix.df$Var2)
colnames(sample.matrix.df) <- c("string", "samples", "value")

int_data <- interactions.robust[ ,c("string", 'type')]

interaction_load <- merge(sample.matrix.df, int_data, by = "string")
interaction_load <- interaction_load %>% group_by(samples, type) %>% summarize(load = sum(value) / length(value))
interaction_load <- reshape2::dcast(samples ~ type, value.var = "load", data = interaction_load)

if (!is.null(sample_data)) {
	interaction_load <- merge(interaction_load, sample_data, by = "samples")
} else {
	interaction_load <- interaction_load
}



if (!is.null(sample_data)) {
	interaction_penetrance <- merge(sample.matrix.df, int_data, by = "string")
	interaction_penetrance <- merge(interaction_penetrance, sample_data, by = "samples")
	interaction_penetrance <- interaction_penetrance %>% group_by(string, type, class) %>% summarize(penetrance = sum(value) / length(value))
	interaction_penetrance <- reshape2::dcast(string + type ~ class, value.var = "penetrance", data = interaction_penetrance)
} else {
	interaction_penetrance <- merge(sample.matrix.df, int_data, by = "string")
	interaction_penetrance <- interaction_penetrance %>% group_by(string, type) %>% summarize(penetrance = sum(value) / length(value))
	interaction_penetrance <- reshape2::dcast(string + type ~ class, value.var = "penetrance", data = interaction_penetrance)
} 


write.table(interaction_load, file = output_load, sep = "\t", quote = F, row.names = F)
write.table(interaction_penetrance, file = output_penetrance, sep = "\t", quote = F, row.names = F)

