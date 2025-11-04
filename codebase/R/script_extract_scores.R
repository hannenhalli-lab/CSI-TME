library(argparse)
parser <- ArgumentParser(description = "A script extract scores from ICA models")
parser$add_argument("--model_file", help = "a .Rds file containing the results of ICA", required = T)
parser$add_argument("--output_file", help = "names of the output file", required = F, default = NULL)


args <- parser$parse_args()

model_file <- args$model_file
ica_sign_correction <- args$ica_sign_correction
output_file <- args$output_file

if (is.null(output_file)) {
	basename_input <- tools::file_path_sans_ext(basename(model_file))
	output_file <- paste("scores_", basename_input, "-Nonaligned.txt", sep = "")
}


library(JADE)

temp_env <- new.env()
load(model_file, envir = temp_env)
#ica_model <- get(ls(temp_env))
ica_model <- temp_env[[ls(temp_env)]]
rm(temp_env)

componentMatix <- ica_model$A
colnames(componentMatix) <- paste("IC.", 1:ncol(componentMatix), sep = "")

#https://stackoverflow.com/questions/2478352/write-table-writes-unwanted-leading-empty-column-to-header-when-has-rownames
#componentMatix <- data.frame("H" = rownames(componentMatix), componentMatix)
sample_names <- gsub("\\.", "-", rownames(componentMatix))
sample_names <- make.names(TCGAutils::TCGAbarcode(sample_names, participant = TRUE), unique = T)
sample_names <- gsub("\\.", "-", sample_names)
componentMatix <- data.frame("H" = sample_names, componentMatix)

message(paste("writing the sample scores for dataset in file:", output_file, sep = " "))
write.table(componentMatix, file = output_file, quote = F, row.names = F, sep = "\t")
