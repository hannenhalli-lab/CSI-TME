library(glue); 
library(dplyr); library(data.table);
library(maftools);  library(magrittr)



selCols <- c("Hugo_Symbol","Entrez_Gene_Id","Center","NCBI_Build","Chromosome","Start_Position",
			"End_Position","Strand","Variant_Classification","Variant_Type","Reference_Allele",
			"Tumor_Seq_Allele1","Tumor_Seq_Allele2","Tumor_Sample_Barcode")
			
#cancerTypes <- c("LAML", "LUAD", "LUSC", "PAAD", "STAD", "SARC", "BLCA","HNSC", "LGG", "LIHC", "PRAD", "THCA")
#cancerTypes <- "LIHC"
cancerTypes <- "glioma.mutation"
for (cancerType in cancerTypes) {
	allFiles <- list.files(glue("GDC/{cancerType}"), recursive = T)
	mafFiles <- grep("maf.gz$", allFiles, value = T)
	mafFiles <- glue("GDC/{cancerType}/{mafFiles}")
	mergedMafFiles <- merge_mafs(mafFiles)
	#dir.create(file.path(glue("MAFs/{cancerType}/")))
	
	#nsynData <- mergedMafFiles@data
	#nsynData <- nsynData %>% subset(., select = selCols)
	
	#synData <- mergedMafFiles@maf.silent
	#synData <- synData %>% subset(., select = selCols)
	
	#mergedMafFiles@data <- nsynData
	#mergedMafFiles@maf.silent <- synData
	
	write.mafSummary(maf = mergedMafFiles, basename = glue("{cancerType}"))
}
