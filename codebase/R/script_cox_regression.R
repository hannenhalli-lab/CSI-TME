library(argparse)
parser <- ArgumentParser(description = "A script to fit cox proportional hazards model")
parser$add_argument("--clinical_file", help = "a tab delimited file containing the clinical data", required = T)
parser$add_argument("--covariates", help = "a comma separated list of co-variates (must be present in clinical_file)", required = T)
parser$add_argument("--strata", help = "at most one variable to be used as strata (must be present in clinical_file)", default = NULL)
parser$add_argument("--feature_file", help = "a tab delimited file containing the predictors (feature by sample)", required = T)
parser$add_argument("--bin_feature", help = "optionally, how many bins the input feature should be binned into (deafult is no binning)",  type = "integer", default = NULL)
parser$add_argument("--samples_along", help = "a character string defining if the feature file contain samples in 'rows' or 'columns'", required = T)
parser$add_argument("--output_file", help = "name of the file to store the output", required = F, default = NULL)
parser$add_argument("--identifier", help = "a name to add for data group", required = F, default = NULL)

args <- parser$parse_args()


suppressPackageStartupMessages( 
	{
		library(dplyr)
		source('codebase/R/functionsForInteractionAnalysis.R')
		library(magrittr)
		library(survival)
		library(survminer)
	} 
	)
	
clinical_file <- args$clinical_file
feature_file <- args$feature_file
samples_along <- args$samples_along
covariates <- args$covariates
strata <- args$strata
output_file <- args$output_file
identifier <- args$identifier
bin_feature <- args$bin_feature

#clinical_file <- "/Users/singha30/DataAndAnalysis/TMEvMutation/TCGA-BRCA_ORACLE/input/TCGA-BRCA_clinical_pam50.txt"
#feature_file <- "/Users/singha30/DataAndAnalysis/TMEvMutation/TCGA-BRCA_ORACLE/input/celltype_ICA/Scores_10/TCGA_BRCA_primary_TPMs_BayesPrism_CAFs_filtered_1.00_ICA_10-aligned.txt"
#samples_along <- "rows"
#covariates <- "age"
#strata <- "stage_class"
#output_file <- NULL
basename_input <- tools::file_path_sans_ext(basename(feature_file))

if (is.null(output_file)) {
	output_file <- paste("survival_cox_coefficients_", basename_input, ".txt", sep = "")
}

covariates <- strsplit(covariates, split = ",") %>% unlist()
covariates_str <- paste(covariates, collapse = " + ")
if (is.null(strata)) {
	formula_str <- paste("Surv(time, status)", "~", "feature", "+", covariates_str, sep = "")
} else {
	#message("***strata is defined")
	if (length(strata) == 1) {
		formula_str <- paste("Surv(time, status)", "~", "feature", "+", covariates_str, "+", "strata(", strata, ")", sep = "")
	} else {
		stop("Error: only a single variable stratification is supported")
	}
}
model_formula <- as.formula(formula_str)
message(model_formula)


#message("reading clinical data****")
clinical_data <- read.table(clinical_file, sep = "\t", header = T)
#message("reading input features data****")
input_data <- read.table(feature_file, sep = "\t", header = T, row.names = 1)


clinical_data <- clinical_data[ ,c('samples', 'time', 'status', covariates, strata)] %>% na.omit
#head(clinical_data)

for (covariate in covariates) {
	if (is.character(clinical_data[[covariate]])) {
		clinical_data[[covariate]] <- factor(clinical_data[[covariate]])
	}
}

if (!is.null(strata)) {
	if (is.character(clinical_data[[strata]])) {
		clinical_data[[strata]] <- factor(clinical_data[[strata]])
		for (covariate in covariates) {
			if (is.numeric(clinical_data[[covariate]])) {
				clinical_data[[covariate]] <- scale(clinical_data[[covariate]], center = T)
			}
		}
	}
} 


if (samples_along == 'columns') {
	input_data <- t(input_data) %>% data.frame()
} else {
	input_data <- input_data %>% data.frame()
}

#message("performing survival analysis")
#head(clinical_data)
survivalAnalysis.list <- list()

for (feature in colnames(input_data)) {
	#feature_data <- input_data[[feature]]		
	feature_data <- input_data[[feature]] %>% set_names(rownames(input_data))
	
	expression.quantiles <- quantile(na.omit(feature_data)) %>% data.frame() %>% t()
	
	feature_data <- data.frame("samples" = names(feature_data), "feature" = feature_data)
	feature_data <- merge(clinical_data, feature_data, by = "samples")
	
	
	if (!is.null(bin_feature)) {
		feature_data$original <- feature_data$feature
		bins <- get.bin.intervals(feature_data$feature, nbins = bin_feature)
		feature_data$feature <- bins
	}
	if (!exists("validate_data")) {
		write.table(feature_data, file = "cox_clinical_data_validation.txt", sep = "\t", row.names = F, quote = F)
		validate_data <- 1
	}
	
	#message(paste("final data has", nrow(feature_data), "rows"))
    #feature_data$age <- scale(feature_data$age, center = T)
	
    HR <- NA; p.val = NA; out <- NA
    tryCatch({
  	  #cox_fit = coxph(Surv(time, status) ~ feature + age + sex, data = feature_data)
	  cox_fit <- coxph(model_formula, data = feature_data)
  	  HR <- summary(cox_fit)$coefficients['feature','coef']
  	  p.val <- summary(cox_fit)$coefficients['feature','Pr(>|z|)']
  	  out <- cbind(feature = feature, hazard = HR, p.value = p.val, expression.quantiles) %>% data.frame()
    },error=function(x){})
		
	if (!is.null(strata)) {
		levels_strata <- unique(feature_data[[strata]])
		cox_fit <- coxph(model_formula, data = feature_data[feature_data[[strata]] == levels_strata[1], ])
    	HR_1 <- summary(cox_fit)$coefficients['feature','coef']
		
		cox_fit <- coxph(model_formula, data = feature_data[feature_data[[strata]] == levels_strata[2], ])
    	HR_2 <- summary(cox_fit)$coefficients['feature','coef']
		
		HR_stratas <- cbind(HR_1, HR_2)
		colnames(HR_stratas) <- levels_strata[1:2]		
		out <- cbind(out, HR_stratas)
	}
	survivalAnalysis.list[[feature]] <- out
}

survivalAnalysis.Df <- do.call("rbind", survivalAnalysis.list) %>% data.frame() %>% na.omit()
survivalAnalysis.Df$p.value <- as.numeric(survivalAnalysis.Df$p.value)
survivalAnalysis.Df$FDR <- p.adjust(survivalAnalysis.Df$p.value, method = "fdr")
survivalAnalysis.Df$hazard <- as.numeric(survivalAnalysis.Df$hazard)
survivalAnalysis.Df$effect <- ifelse(survivalAnalysis.Df$hazard > 0, "Pro-Tumor", "Anti-Tumor")

if(!is.null(identifier)) {
	survivalAnalysis.Df$group <- identifier
}
write.table(survivalAnalysis.Df, output_file, sep = "\t", row.names = F, quote = F)
message(paste("done!. cox output stored in", output_file))


#surv_object <- Surv(clinical_data$time, clinical_data$status)
#fit <- survfit(surv_object ~ feature, data = clinical_data)
#p <- ggsurvplot(fit)

