args <- commandArgs(trailingOnly = TRUE)

######Validation######
#project <- "TCGA-LUSC"; target = "NSCLC_Sorafenib_2"

project <- args[1]
target <- args[2]

for (target in targets) {
	projections_input <-  paste("input/", project, "/projections", sep = "")
	projections_output <-  paste("output/", project, "/projections/", target, sep = "")

	dir.create(projections_output, recursive = T)
	dir.create(projections_input, recursive = T)
}
