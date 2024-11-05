library(TCGAbiolinks); library(glue); library(dplyr); library(maftools)

#cancer <- "CHOL"

c("LAML","ACC","BLCA","LGG","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC",
	"KICH","KIRC","KIRP","LIHC","LUAD","LUSC","DLBC","MESO","OV",
	"PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THYM","THCA","UCS","UCEC","UVM")
	
for (cancer in cancerTypes) {
	query <- GDCquery(
	    project = glue("TCGA-{cancer}"),
	    data.category = "Simple Nucleotide Variation", 
	    access = "open",
	    data.type = "Masked Somatic Mutation", 
	    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
	)
	GDCdownload(query = query)
	maf.data <- GDCprepare(query = query, summarizedExperiment = F)
	write.table(maf.data, file = glue("TCGAdata/{cancer}.maf"), sep = "\t", quote = F, row.names = F)
}








########################obsolete########################

#cna.data <- GDCprepare(query = query)
#cna.data <- SummarizedExperiment::assays(cna.data)$copy_number
#colnames(cna.data) <- gsub(",.+", "", colnames(cna.data))