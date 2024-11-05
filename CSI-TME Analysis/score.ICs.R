library(Seurat); library(glue); library(dplyr)


cohort <- "tcga"; idh = "mut"
glioma <- readRDS('all_merged_after_markers_n_with_anno.RDS')
signatureFile <- glue('../ICA-JADE-10.signatureGenes.{cohort}.IDH{idh}--filtered-noNorm-optimalRank.Rda')

load(signatureFile)



glioma.idh <- subset(x = glioma, subset = IDH == "IDH Mut")
cellTypes <- names(signature.list)


top = "negative"
for (cellType in cellTypes) {
	geneset <- signature.list[[cellType]]
	geneset.new <- list()
	for (dims in names(geneset)) {
		if (top == "positive") {
			l1 <- geneset[[dims]]$positive
		}
		if (top == "negative") {
			l1 <- geneset[[dims]]$negative
		}
		if (length(l1) > 2) {
			geneset.new[[dims]] <- l1
		}
	}
	
	cellType.subset <- cellType
	if (cellType == "Bcell") {
		cellType.subset <- "B cell"
	}
	if (cellType == "Tcell") {
		cellType.subset <- "T cell"
	}
	
	glioma.idh.cellType <- subset(x = glioma.idh, subset = Anno == cellType.subset)
	
	glioma.idh.cellType <- AddModuleScore(
	  object = glioma.idh.cellType,
	  features = geneset.new,
	  name = glue("ICA.{cellType}")
	)
	
	for (index in 1:length(geneset.new)) {
		pdf(file = glue("ICA.component{index}.{top}.{cellType}.{cohort}.IDH{idh}.pdf"))
		f.plot <- FeaturePlot(object = glioma.idh.cellType, features = glue("ICA.{cellType}{index}"), min.cutoff = 0)
		print(f.plot)
		dev.off()
	}
	message(glue("done for {cellType}"))
}



################################################################
################################################################
################################################################

pdf(file = glue("codeletion.{cellType}.{cohort}.IDH{idh}.pdf"))
FeaturePlot(object = glioma.idh.cellType, features = glue("NMF.{cellType}.state{index}"), min.cutoff = 0)
dev.off()

