# Overview

<img width="576" alt="image" src="https://github.com/user-attachments/assets/07e9da7b-343d-43f8-b110-acd7e53c9dcf">

CSI-TME - A computational pipeline to infer cell state interaction network using bulk transcriptomic and clinical data from cancer patients. 
Typically, cell states in tumor microenvironments (TME) are only investigated using single cell datasets. However, there are two key limitations – 

- Single cell datasets are far less abundant as well as have smaller cohort sizes than bulk transcriptomic datasets making it difficult to make reliable clinical predictions.

- Secondly, single cell datasets are rarely linked to the clinical outcomes collected in a principled fashion.

To overcome these limitations, CSI-TME uses a combination of supervised and unsupervised learning approaches to first the estimate cell state composition in the tumor microenvironment, and then model the clinical data based on the joint activity of the pairs of cell states.

For details, please check our preprint @ https://www.biorxiv.org/content/10.1101/2024.10.29.620901v1

# Running CSI-TME

There are two main use cases for CSI-TME

1. Minining prognostic cell state interactions in a new cohort.
2. Using pre-trained model to score your samples for pro- and anti-tumor interactions.

A step-by-step tutorial is provided below for each of these steps 

## Required input
- Bulk gene expression data from a large cohort of cancer patients (Gene x Sample matrix). We recommend atleast 200 samples for sufficient statistical power.

- Matching clinical data with atleast the following columns

<img width="468" height="141" alt="clinical" src="https://github.com/user-attachments/assets/6a639c23-4b11-41cb-a437-13b806fbbe1c" />

status 0 => Alive; 1 => Death

### Step 1. Create directories to control input and output flow
From Unix shell, script_path_training.R to create the directories where input/output data will be stored; for eg.

``` Rscript codebase/R/script_path_training.R Cohort_1 ```

This will create two directories; **input/Cohort_1** and **output/Cohort_1**. 

- If you have already deconcoved the cell type specific gene expression profiles, please store them at **input/Cohort_1/deconvolution** with a spearate file name for each cell type (for eg. input/Cohort_1/deconvolution/Tcell.txt). 

- Otherwise, place the bulk gene expression file  (Gene x Samples) at the **input/Cohort_1/genes/expression_data.txt** file and perform the cell type specific deconvolution using 

``` bash codebase/bash/deconvolve.sh Cohort_1 ```

This will perform deconvolution using a combination of cibersortX and CODEFACS and the deconvoved data will be placed in **input/Cohort_1/deconvolution** directory. We strongly encourage to try a couple of different deconvolution algorithams. 

### Step 2. Detection cell state interactions

After deconvolution, run the following command

 ``` bash codebase/bash/CSI-TME.sh Cohort_1 ```

This will perform create ICA factorization of the cell type specific gene expression profiles and screen for prognostic cell state combinations of the IC pairs using Cox regression, followed by 10 bootstraps to deduce crossvalidation accuracy. The final output will be stored in file called CSI-TME_significant_crossvalidation.txt in the directory named as **input/Cohort_1**. Users are free to set select most significant interactions of their interest using FDR thresholds or crossvalidation accuracy thresholds. 

### Step 3. Validate cell state interactions in an independent cohort

If users have access to the bulk transcriptomic datasets with matched clinical outcomes from more than one large cohorts, then remaining cohorts can be used for the purpose of cross-cohort validation as follows

``` Rscript codebase/R/script_path_validation.R Cohort_2 ``` 

place the bulk gene expression data in **input/Cohort_1/validation/Cohort_2/genes/expression_data.txt** and clinical data in **input/Cohort_1/validation/Cohort_2/clinical_data.txt** for new cohort and run 

``` bash codebase/bash/ica_project_new.sh Cohort_1 Cohort_2 ``` 

This will use the gene weights from ICA models fitted Cohort_1 to project the bulk data from Cohort_2 

Lastly -

``` bash codebase/bash/CSI_validate.sh Cohort_1 Cohort_2 ```

This step will calculcate the hazard ratios using the clinical data for new Cohort_2. The output from this step will be stored in file **output/Cohort_1/validation/Cohort_2/CSI_discovery_validation.txt**. 

### Step 4. Exploratory analysis of cell state interactions

Finally, we provide some options to help interprettation and additional exploratory analysis  as follows


#### Extract marker genes of Independent Components

``` bash codebase/R/Rscript script_extract_markers.R --cohort_name Cohort_1 ```

This will create a file **output/Cohort_1/marker_genes.txt** listing cell type, their ICs, and corresponding marker genes. These marker genes can be used to interpret the CSIs

#### Extract scores of Independent Components across samples in discovery cohort, which correspond to the cell states

``` bash codebase/R/Rscript script_extract_cellstates.R --cohort_name Cohort_1 ```

#### Extract the interaction penetrance and interaction load

* After identifying and validating CSIs, one can extract the interaction penetrance (calculated per interaction), and interaction load (calculated per sample) in the discovery or any new cohort. 

* If there are more than one kinds of samples (for instance, responders and non-responders) then both interaction penetrance and interation load can be calculated separately for each sample group by providing the sample classifications with ``` --sample_file ``` argument. Otherwise, skip using ``` --sample_file ``` argument

``` bash codebase/R/Rscript script_penetrance_load.R --discovery_cohort Cohort_1 --target_cohort Cohort_2 --is_cohort validation --sample_file input/Cohort_1/validation/Cohort_2/sample_class.txt ```

This will generate two files - 
* output/Cohort_1/interaction_load_Cohort_2.txt
* output/Cohort_1/interaction_penetrance_Cohort_2.txt




_If you have any questions, please reach out at arashdeep.singh@nih.gov or sridhar.hannenhalli@nih.gov_


