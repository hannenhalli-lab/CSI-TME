library(argparse)

# Set up argument parser
parser <- ArgumentParser(description = "Perform regression analysis")
parser$add_argument("--cohort_name", type = "character", help = "Name of the discovery cohort", required = TRUE)
parser$add_argument("--project", type = "character", help = "name of the project", required = TRUE)


args <- parser$parse_args()
cohort_name <- args$cohort_name
project <- args$project


library(dplyr); library(magrittr)
source("codebase/R/functionsForInteractionAnalysis.R")


ica_model_files <- list.files(cohort_name, full.name = T)
sample.names <- NULL

for (ica_model_file in ica_model_files) {
	cell_type <- gsub("\\.Rda", "", basename(ica_model_file))
	load(ica_model_file)
  	pcaRes <- ica_model
 	component_matix <- ica_model$A
	gene_matrix <- ica_model$S
	colnames(component_matix) <- colnames(gene_matrix)
	
	component_file <- paste("output/", project, "/cellstates/", cell_type, ".txt", sep = "")
	gene_file <- paste("output/", project, "/genes/", cell_type, ".txt", sep = "")
	
	write.table(component_matix, file = component_file, sep = "\t", quote = F)
	#write.table(gene_matrix, file = gene_file, sep = "\t", quote = F)
}
message("extracted cell states and gene loadings are succesfully written")


