#project <- "TCGA-IDH"
project <- commandArgs(trailingOnly = TRUE)



gene_input <- paste("input/", project, "/genes", sep = "")
cellstates_input <- paste("input/", project, "/cellstates", sep = "")
deconvolution_input <- paste("input/", project, "/deconvolution", sep = "")
projections_input <-  paste("input/", project, "/projections", sep = "")


gene_output <- paste("output/", project, "/genes", sep = "")
cellstates_output <- paste("output/", project, "/cellstates", sep = "")
models_output <- paste("output/", project, "/ica_models", sep = "")
projections_output <-  paste("output/", project, "/projections", sep = "")


dir.create(gene_input, recursive = T)
dir.create(cellstates_input, recursive = T)
dir.create(deconvolution_input, recursive = T)
dir.create(projections_input, recursive = T)


dir.create(gene_output, recursive = T)
dir.create(cellstates_output, recursive = T)
dir.create(models_output, recursive = T)
dir.create(projections_output, recursive = T)








