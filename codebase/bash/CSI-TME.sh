#!/bin/bash
#BRCA,BLCA,ESCA,HNSC,KIRC,LIHC,LUAD,STAD
#project="TCGA-BLCA"
project=$1

input_directory="input/${project}/deconvolution/"
clinical_file="input/${project}/clinical_data.txt"
n_comp=10
covariates="age,sex"
strata="stage_class"

output_ica="output/${project}/ica_models"
output_CSI="output/${project}/CSI-TME_interactions.txt"
output_bootstrapped="output/${project}/CSI-TME_significant_crossvalidation.txt"


echo "finding ICA decompositions..."
Rscript codebase/R/script_ica.R --input_directory $input_directory --n_comp $n_comp --output_dir $output_ica

echo "performing pairwise comparison of cell types..."

if [[ -z "$strata" ]]; then
	Rscript codebase/R/CSI-TME_main.R \
		--ica_model_directory $output_ica \
		--clinical_file $clinical_file \
		--covariates $covariates \
		--output_file $output_CSI \
		--is_tcga
else
	Rscript codebase/R/CSI-TME_main.R \
		--ica_model_directory $output_ica \
		--clinical_file $clinical_file \
		--covariates $covariates \
		--strata $strata \
		--output_file $output_CSI \
		--is_tcga
fi


echo "performing crossvalidation using bootstrapping..."

if [[ -z "$strata" ]]; then
	Rscript codebase/R/script_bootstrapping.R \
		--ica_model_directory $output_ica \
		--clinical_file $clinical_file \
		--interaction_file $output_CSI \
		--covariates $covariates \
		--output_file $output_bootstrapped \
		--is_tcga
else
	Rscript codebase/R/script_bootstrapping.R \
		--ica_model_directory $output_ica \
		--clinical_file $clinical_file \
		--interaction_file $output_CSI \
		--covariates $covariates \
		--strata $strata \
		--output_file $output_bootstrapped \
		--is_tcga
fi

