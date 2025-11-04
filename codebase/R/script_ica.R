
library(argparse)

# Set up argument parser
parser <- ArgumentParser(description = "Perform regression analysis")
parser$add_argument("--input_directory", type = "character", help = "a directory with cell type specific gene expression file (genes by sample) ", required = TRUE)
parser$add_argument("--n_comp", type = "numeric", help = "how many ICs to extract?", required = F, default = 10)
parser$add_argument("--output_dir", type = "character", help = "a directory to save the ICA models", required = T) #####default doesnt store the fitted models

args <- parser$parse_args()
input_directory <- args$input_directory
n_comp <- args$n_comp
output_dir <- args$output_dir

message(input_directory)

input_files <- list.files(input_directory, full.name = T)



library(MineICA); library(dplyr)
library(magrittr); library(glue)
library(doParallel); library(doSNOW)
source('codebase/R/functionsForInteractionAnalysis.R')


method = "JADE"; max.dim = 10; dims = n_comp
###progress monitor solution from  https://stackoverflow.com/questions/5423760/how-do-you-create-a-progress-bar-when-using-the-foreach-function-in-r


#dims.cellTypes <- c(Bcell = 42, Endothelial = 40, Malignant = 59, Myeloid = 47, Oligos = 52, Stromal = 44, Tcell = 45)
set.seed(123456)
ica.res.list <- list()
for (input_file in input_files) {
	extension <- tools::file_ext(input_file)
	if (any(extension %in% c("txt", "tsv"))) {
		inputFile <- read.table(input_file, sep = "\t", header = T)
	}
	if (any(extension %in% c("csv"))) {
		inputFile <- read.csv(input_file,  header = T)
	}
	cell_type <- gsub("\\.txt", "", basename(input_file))
	expression.data <- preProcessInput(inputFile, removeNA = T, removeSD = T, removeMedian = T, logT = T, normalize = F, zscore = T, positive = F)
	if (dims == "MSTD") {
		max.dim <- dims.cellTypes[[cell_type]]
	}
	if (method == "JADE") {
		ica_model <- runICA(method = "JADE", X = expression.data, nbComp = max.dim, alg.type = "parallel", maxit = 20000, tol = 10^-6)
			
	} else {
		###this is slow due to multiple runs for each decomposition
		cl <- makeCluster(future::availableCores() - 4)
		registerDoSNOW(cl)
		ica_model <- clusterFastICARuns.modified(X = expression.data, nbComp = max.dim, alg.type = "deflation", nbIt = 50, funClus = "hclust", method = "average")
		stopCluster(cl);
	}
	ica.res.list[[cell_type]] <- ica_model
	print(glue("done for {cell_type}"))
	output_file = paste(output_dir, "/", cell_type, ".Rda", sep = "")
	save(ica_model, file = output_file)
}
message("ICA models are fitted and saved successfully...")

#save(ica.res.list, file = output_file)


