library(argparse)

# Set up argument parser
parser <- ArgumentParser(description = "Perform regression analysis")
parser$add_argument("--ica_model_directory", type = "character", help = "a directory with .Rda files of the fitted ICA models", required = TRUE)
parser$add_argument("--new_expression_file", type = "character", help = "a file containing the clinical data", required = TRUE)
parser$add_argument("--target_cell_type", type = "character", help = "deisred cell type name for ica model to use, by default, all cell types in model directory are used", required = F, default = NULL)
parser$add_argument("--zscore_output", action = "store_true", help = "should the IC scores be z-scored relative to a null model", required = F, default = F)
parser$add_argument("--output_dir", type = "character", help = "a directory to the projected ICA models", required = T) #####default doesnt store the fitted models

suppressPackageStartupMessages( 
	{
		library(dplyr); library(magrittr); 
		library(preprocessCore); library(glue); 
		source('codebase/R/functionsForInteractionAnalysis.R')
	}
	)




project_ica <- function(new_data, ica_model) {
	S <- ica_model$S
	common_genes <- intersect(rownames(S), rownames(new_data))
	
	new_data <- new_data[common_genes, ]
	discovery_rotation <- S[common_genes, ]
	warning(paste("ICA projection of new data is based on ", length(common_genes), "common genes"))
	projection <- t(new_data) %*% discovery_rotation
    projection <- apply(projection, 2, function(x) scales::rescale(x, to = c(-1, 1)))
	ica <- list()
	ica$A <- projection
	ica$S <- NA
	ica$W <- NA
	ica$Xmu <- NA
	#return(projection)
	return(ica)
}


progress_monitor <- function(n_iter, current_iter, percent_offset = 5) {
	percent_complete <- (current_iter / n_iter) * 100
  
  	last_printed <- floor((current_iter - 1) / (n_iter * percent_offset / 100)) * percent_offset
	next_print <- floor(percent_complete / percent_offset) * percent_offset
  
	if (next_print > last_printed) {
		cat(sprintf("Progress: %d%% completed\n", next_print))
		flush.console()
	}
}


args <- parser$parse_args()
ica_model_directory <- args$ica_model_directory
new_expression_file <- args$new_expression_file
target_cell_type <- args$target_cell_type
output_dir <- args$output_dir
zscore_output <- args$zscore_output


input_data <- read.table(new_expression_file, sep = "\t", header = T, row.names = 1)
ica_models_files <- list.files(ica_model_directory, full.name = T)

if (!is.null(target_cell_type)) {
	ica_model_files = grep(target_cell_type, ica_models_files, value = T)
}

input_data <- preProcessInput(input_data, removeNA = T, removeSD = T, removeMedian = T, logT = F, normalize = T, zscore = T, positive = F)

ica.res.list <- list()
set.seed(12345); n_iterations = 100
for (model_file in ica_models_files) {
	cell_type <- gsub("\\.Rda", "", basename(model_file))
	message(model_file)
	
	temp_env <- new.env()
	load(model_file, envir = temp_env)
	#ica_model <- get(ls(temp_env))
	ica_model <- temp_env[[ls(temp_env)]]
	rm(temp_env)
	
	
	
	#A <- ica_model$A
	projected_model <- project_ica(new_data = input_data, ica_model = ica_model)
	
	#output_file <- paste(output_dir, "/scores_", "celltype", ".txt", sep = "")
	#write.table(projected_model$A, file = output_file, sep = "\t", quote = F)
	
	output_file <- paste(output_dir, "/ica_models/", cell_type, ".Rda", sep = "")
	save(projected_model, file = output_file)
	
	
	
	output_file <- paste(output_dir, "/ica_scores/", cell_type, ".txt", sep = "")
	A <- projected_model$A
	if (zscore_output) {
		null_space <- list()
		message(paste("generating randomized ICA model for ", cell_type))
		for (i in 1:n_iterations) {
			random_vector <- sample(rownames(input_data))
			random_data <- input_data
			rownames(random_data) <- random_vector
			null_model <- project_ica(new_data = random_data, ica_model = ica_model)
			null_space[[i]] <- null_model$A
			#progress_monitor(n_iter = n_iterations, current_iter = i, percent_offset = 5)
		}
	
		#null_mean_df <- lapply(null_space, function(x) apply(x, 2, mean)) %>% do.call("rbind", .)
		#null_mean <- apply(null_mean_df, 2, mean)
	
		#null_sd_df <- lapply(null_space, function(x) apply(x, 2, sd)) %>% do.call("rbind", .)
		#null_sd <- apply(null_sd_df, 2, mean)
	
		null_comps <- do.call("rbind", null_space)
		null_mean <- apply(null_comps, 2, mean, na.rm = T)
		null_sd <- apply(null_comps, 2, sd, na.rm = T)
	
		zscored_projections <- sapply(1:length(null_mean), function(index) {(A[ ,index] - null_mean[index]) / null_sd[index]})
		colnames(zscored_projections) <- names(null_mean)
		write.table(zscored_projections, file = output_file, sep = "\t", quote = F)
		message(paste("done for ", cell_type))
	} else {
		write.table(A, file = output_file, sep = "\t", quote = F)
	}
	#output_file <- paste(output_dir, "/scores_", "celltype", ".txt", sep = "")
	#write.table(projection$A, file = output_file, sep = "\t", quote = F)
}
