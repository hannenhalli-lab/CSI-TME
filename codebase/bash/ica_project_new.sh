#!/bin/bash
#!/bin/bash

project=$1; 
target=$2;

ica_models=output/${project}/ica_models

new_expression_file=input/${project}/projections/${target}_expression.txt

projected_models_dir=output/${project}/projections/${target}/ica_models
projected_scores_dir=output/${project}/projections/${target}/ica_scores ###these are z-scored relative to null space

mkdir -p ${projected_models_dir}
mkdir -p ${projected_scores_dir}


echo "projecting data onto the ICA space of ${project}"
Rscript codebase/R/ica_script_project_new.R \
	--ica_model_directory $ica_models \
	--new_expression_file $new_expression_file \
	--output_dir "output/${project}/projections/${target}"

