#!/bin/bash

####requirements, input directory where deconvovled files are located.
project_name="TCGA-COAD"
base_dir="./"

expression_file="${base_dir}/input/${project_name}/genes/${project_name}_primary_TPMs.txt"
clinical_file="${base_dir}/input/${project_name}/clinical_data.txt"

output_dir="${base_dir}/output/${project_name}/genes/"


covaraites="age"
strata="stage_class"


mkdir -p ${output_dir}

output_file=$(basename $expression_file)
output_file="genes_prognosis_${output_file}"
output_file=${output_dir}${output_file}


#####pre processing of the input file to normalize, log transform and then bring in the desired sample (rows) by genes (columns) format
echo "bash --> prcesssing bulk gene expression file for ORACLE analysis.."

processed_file=$(mktemp)
Rscript R/script_process_geneExpression.R  \
	--input_file $expression_file \
	--log_transform yes \
	--quantile_normalize yes \
	--remove_median 1 \
	--output_file $processed_file \
	--show_message yes


header_file=header_file_${project_name}

: '
####this code was used if when the input expression file had one less column in the header (typicall for R generated tab files)
####but using script_process_geneExpression.R already fixes this issue
echo "bash --> processing..."
touch $header_file
echo "H" > $header_file
csvtk headers $processed_file | tr "-" "\." | tr [:space:] "\n" >> $header_file


tail -n +2 $processed_file | csvtk -H tab2csv > ${processed_file}_noHeader 
csvtk -H transpose ${processed_file}_noHeader  > ${processed_file}_noHeader_tranposed ####
paste -d ',' header_file.txt  ${processed_file}_noHeader_tranposed | csvtk csv2tab >  ${processed_file}_noHeader_tranposed_headed
'

temp_file=$(mktemp)

csvtk -t transpose $processed_file > ${processed_file}_transposed

bash/add_tcga_samples.sh ${processed_file}_transposed H | csvtk -t cut -f -H > $temp_file

if [ -n "$strata" ]; then
	Rscript R/script_cox_regression.R --clinical_file $clinical_file --covariates $covaraites --bin_feature 3 --strata $strata --feature_file $temp_file --samples_along rows --output_file $output_file --identifier bulk
else
	Rscript R/script_cox_regression.R --clinical_file $clinical_file --covariates $covaraites --bin_feature 3 --feature_file $temp_file --samples_along rows --output_file $output_file --identifier bulk
fi


