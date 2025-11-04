library(argparse)

# Set up argument parser
parser <- ArgumentParser(description = "extract marker genes of the cell states (ICs)")
parser$add_argument("--cohort_name", type = "character", help = "name of the discovery cohort", required = TRUE)

args <- parser$parse_args()
cohort_name <- args$cohort_name

output_connection <- paste("output/", cohort_name, "/IC_marker_genes.txt")

library(dplyr); library(magrittr);
library(glue);

get_genes <- function(index, S) {
	x <- S[ ,index]
	pos <- x[x > mean(x) + list.cutoffs[index]] %>% names()
	neg <- x[x < mean(x) - list.cutoffs[index]] %>% names()
	
	pos <- data.frame(Gene = pos, direction = "positive")
	neg <- data.frame(Gene = neg, direction = "negative")
	
	#genes <- list("positive" = pos, "negative" = neg)
	genes <- rbind(pos, neg)
	return(genes)
}

genesets_ica <- function(ica_model, cutoff = 2.5) {
	S <- ica_model$S
	stdev <- apply(S, 2, sd)
	
	list.cutoffs <- cutoff * stdev
	
	gene_sets <- sapply(1:length(list.cutoffs), function(index) get_genes(index, S), simplify = F)
	gene_sets <- lapply(gene_sets, function(x) data.frame(Gene = x))
	return(gene_sets)
}
#####select only the

model_directory <- paste("input/", cohort_name, "/ica_models/")

all_files <- list.files(model_directory, full.names = T)

signature_list <- list()
for (model_file in ica_models_files) {
	cell_type <- gsub("\\.Rda", "", basename(model_file))
	message(model_file)
	
	temp_env <- new.env()
	load(model_file, envir = temp_env)
	#ica_model <- get(ls(temp_env))
	ica_model <- temp_env[[ls(temp_env)]]
	rm(temp_env)
	
	gene_list <- genesets_ica(ica_model, cutoff = 2.5) %>% set_names(., paste("IC", 1:length(.), sep = "."))
	
	marker_genes <- bind_rows(gene_list, .id = "IC")
	signature_list[[cellType]] <- gene_list
	
}

signature_df <- bind_rows(signature_list, .id = "cell_type")


message(paste("writing output to", output_connection))
write.table(signature_df, file = output_connection, sep = "\t", quote = F, row.names = F)


