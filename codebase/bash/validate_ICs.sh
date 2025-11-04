#!/bin/bash

project="TCGA-IDH";
validation="CGGA"

new_input_dir=input/${project}/validation/$validation/ ####a path where deconvovled data for validation cohort is placed
output_random_discovery="output/${project}/randomized/ica_scores" ####a path to store the IC scores of randomized discovery cohort
output_random_projections="output/${project}/randomized/projections/$validation" ####a path to store the projections of validation cohort using random IC models
mkdir -p $output_random_projections

max_iterations=10

for ((i=1; i<11; i++)); do 
	ica_models_dir=output/${project}/randomized/ica_models/iteration${i}
	
	for ica_model in ${ica_models_dir}/*; do
		echo $ica_model
		cell_type=$(echo $ica_model | sed 's/.*\///' | sed 's/.Rda//')
		new_expression_file=${new_input_dir}/deconvolution/${cell_type}.txt
		
		mkdir -p ${output_random_projections}/iteration${i}/ica_models
		mkdir -p ${output_random_projections}/iteration${i}/ica_scores
		
		mkdir -p ${output_random_discovery}/iteration${i}/
		
		echo "extracting cell state scores of the randomized data"
		Rscript codebase/R/script_extract_scores.R \
			--model_file $ica_model \
			--output_file ${output_random_discovery}/iteration${i}/${cell_type}.txt
				
		echo "projecting data onto the ICA space of ${project}"
		Rscript codebase/R/ica_script_project_new.R \
			--ica_model_directory ${ica_models_dir} \
			--new_expression_file $new_expression_file \
			--target_cell_type $cell_type \
			--output_random_projections "${output_random_projections}/iteration${i}" \
		
		echo "created random projections for $cell_type & iteration $i"
	done
done


#####perform survival analysis of the randmized discovery cohort data######
for ((i=1; i<11; i++)); do 
	score_directory=${output_random_discovery}/iteration${i}/
	survival_directory=${output_random_discovery}/../prognosis/iteration${i}/
	mkdir -p $survival_directory
	for score_file in $score_directory/*; do 
		cell_type=$(echo $score_file | sed 's/.*\///' | sed 's/.txt//')
		Rscript codebase/R/script_cox_regression.R \
			--clinical_file $new_input_dir/../../clinical_data.txt \
			--covariates "age,sex" \
			--feature_file $score_file \
			--samples_along 'rows' \
			--output_file ${survival_directory}/${cell_type}.txt \
			--identifier "random${i}"
	done
done

: '

#####perform survival analysis of the projected validation cohort data######
for ((i=1; i<11; i++)); do 
	score_directory=${output_random_projections}/iteration${i}/ica_scores
	mkdir ${output_random_projections}/iteration${i}/prognosis
	for score_file in $score_directory/*; do 
		cell_type=$(echo $score_file | sed 's/.*\///' | sed 's/.txt//')
		Rscript codebase/R/script_cox_regression.R \
			--clinical_file $new_input_dir/clinical_data.txt \
			--covariates "age,sex" \
			--feature_file $score_file \
			--samples_along 'rows' \
			--output_file ${output_random_projections}/iteration${i}/prognosis/${cell_type}.txt \
			--identifier "random${i}"
	done
done
'



