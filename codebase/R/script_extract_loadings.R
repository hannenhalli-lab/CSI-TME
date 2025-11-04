library(argparse)
parser <- ArgumentParser(description = "A script extract scores from PCA/ICA/NMF")
parser$add_argument("--factorization_file", help = "a .Rds file containing the results of ICA/NMF/PCA", required = T)
parser$add_argument("--factorization_type", help = "whether input file is ICA, NMF, or PCA", required = T)
parser$add_argument("--ica_sign_correction", action = "store_true", help = "if factorization_type is ICA, should the components be aligned positively with maximum kurtosis", required = F)
parser$add_argument("--output_file", help = "names of the output file", required = F, default = NULL)


args <- parser$parse_args()

factorization_file <- args$factorization_file
factorization_type <- args$factorization_type
ica_sign_correction <- args$ica_sign_correction
output_file <- args$output_file

if (is.null(output_file)) {
	if (ica_sign_correction) {
		basename_input <- tools::file_path_sans_ext(basename(factorization_file))
		output_file <- paste("scores_", basename_input, "-aligned.txt", sep = "")
	} else {
			basename_input <- tools::file_path_sans_ext(basename(factorization_file))
			output_file <- paste("scores_", basename_input, "-Nonaligned.txt", sep = "")
		}
}

scores_list <- list()
if (factorization_type == "PCA") {
	library(FactoMineR)
	
	factorization <- readRDS(factorization_file)
	loadingMatix <- factorization$var$coord
}

if (factorization_type == "NMF") {
	library(NMF)
	
	factorization <- readRDS(factorization_file)
	#componentMatix <- factorization@fit@H %>% t()
}

if (factorization_type == "ICA") {
	library(JADE)
	
	factorization <- readRDS(factorization_file)
	
	componentMatix <- factorization$S
	loadingMatix <- factorization$A
	
	if(ica_sign_correction) {
		library(parameters)
		
		k_values <- apply(componentMatix, 2, function(x) parameters::skewness(x)[1])
		k_values <- unlist(k_values)
		k_signs <- sign(k_values)
		
		#k_signs <- apply(componentMatix, 2, function(x) sign(parameters::skewness(x)[1]))
		
		for (i in 1:length(k_signs)) {
			####print the skewness values to get a sense
			cat(sprintf("skewness of %s is %f \n", paste("IC",i, sep = "."), k_values[i]))
			flush.console()
			
			k_sign <- k_signs[i]
			if (k_sign != 0 ){
				loadingMatix[ ,i] <- loadingMatix[ ,i] * k_sign
			}
		}
	} 
		
	rownames(loadingMatix) <- names(factorization$Xmu)
	colnames(loadingMatix) <- paste("IC", 1:ncol(loadingMatix), sep = ".")
}

message(paste("writing the loadings scores for dataset in file:", output_file, sep = " "))

loadingMatix <- data.frame("H" = rownames(loadingMatix), loadingMatix) 
#https://stackoverflow.com/questions/2478352/write-table-writes-unwanted-leading-empty-column-to-header-when-has-rownames
write.table(loadingMatix, file = output_file, quote = F, row.names = F, sep = "\t")
