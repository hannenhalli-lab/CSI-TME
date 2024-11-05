library(Seurat); library(glue); 
library(tibble); library(dplyr)
library(amap); library(usedist)

library(doParallel);



###important functions to maniputate distances######
####https://cran.r-project.org/web/packages/usedist/readme/README.html
#### used for current purposes 


cohort <- "tcga"; idh = "mut"
seuratFile <- "/data/singha30/Deconvolution/codefacs/newDeconv/Seurat/seurat.glioma.idhMut.rds"
signatureFile <- glue('/data/singha30/Deconvolution/codefacs/newDeconv/ICA-JADE-10.signatureGenes.{cohort}.IDH{idh}--filtered-noNorm-optimalRank.Rda')


glioma.idh <- readRDS(seuratFile)
load(signatureFile)
message("done loading the seurat file***\n")


myargs = commandArgs(trailingOnly=TRUE)
direction = myargs[1]
outputDir <- "/data/singha30/Deconvolution/codefacs/newDeconv"



#####compute distances matrix using Dist function with parallel option
message("now computing distance matrices***\n")

size.parallel <- future::availableCores()
dist_matrix_list <- list()
cellTypes <- names(signature.list)

for (cellType in cellTypes) {
	cellType.subset <- cellType
	if (cellType == "Bcell") {
		cellType.subset <- "B cell"
	}
	if (cellType == "Tcell") {
		cellType.subset <- "T cell"
	}
	
	glioma.idh.cellType <- subset(x = glioma.idh, subset = Anno == cellType.subset)
    pca_coords <- Embeddings(glioma.idh.cellType, "pca")
	pca_coords <- pca_coords[ ,1:10]
	dist_matrix <- Dist(pca_coords, method = "euclidean", nbproc = size.parallel, diag = T, upper = T)
	dist_matrix <- as.matrix(dist_matrix)
	dist_matrix_list[[cellType]] <- dist_matrix
	message(glue("computed distance matrix for {cellType}"))
}
#table(glioma.idh@meta.data$Anno)

#######first, add the signature scores

distances_celltype_list <- list()
#direction = "negative"
for (cellType in cellTypes) {
	message(glue("now scoring signatures for {cellType}***\n"))
	geneset <- signature.list[[cellType]]
	geneset.new <- list()
	for (dims in names(geneset)) {
		if (direction == "positive") {
			l1 <- geneset[[dims]]$positive
		}
		if (direction == "negative") {
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
	
	suppressWarnings({
		glioma.idh.cellType <- AddModuleScore(
		  object = glioma.idh.cellType,
		  features = geneset.new,
		  name = names(geneset.new)
		)
	})
	
		####the following two steps are needed becuase AddModuleScore score adds extra numerics at the end of provided module names####
	module.names.indices <- grep(glue("Dim.\\w+"), colnames(glioma.idh.cellType@meta.data))
	colnames(glioma.idh.cellType@meta.data)[module.names.indices] <- names(geneset.new)
	
	
	meta.data <- glioma.idh.cellType@meta.data[ ,module.names.indices]
	meta.data <- apply(meta.data, 2, function(x) ifelse(x > 0, "Expressed", "Remaining")) %>% data.frame()
	
	####too many malignant cells make the following calculations difficult. Subsampling used for the purpose of efficiency.
	if (cellType == "Malignant" | cellType == "Myeloid") {
		reduce.fold <- 5
		n <- nrow(meta.data)
		sample.size <- floor(n/reduce.fold)
		percent.retained <- round((1 / reduce.fold) * 100, 2)
		indices <- sample(1:n, sample.size, replace = F)
		meta.data <- meta.data[indices, ]
		message(glue("retained {percent.retained}% of the cells for {cellType}***\n"))
	}	
	dist_matrix <- dist_matrix_list[[cellType]]
	distances_list <- list()
	
	message(glue("now extracting distances ICs of {cellType}***\n"))
	
	#size.loop <- length(names(meta.data))
	
	#cl <- makeCluster(future::availableCores())
	#registerDoParallel(cl)
	#registerDoMC(cores=future::availableCores())
	
	#registerDoSNOW(cl)
	#pb <- txtProgressBar(max = size.loop, style = 3)
	#progress <- function(n) setTxtProgressBar(pb, n)
	#opts <- list(progress = progress)
	
	
    #distances_list <- foreach(index = 1:size.loop, .packages=c('usedist', 'dplyr', 'magrittr', 'tibble' ,'glue')) %dopar% {
	for (index in 1:length(names(meta.data))){
		sig <- names(meta.data)[index]
	    membership <- meta.data %>% select(., select = paste(sig)) %>% rownames_to_column(., "cells")
		membership <- split(membership, membership$select) %>% lapply(., function(x) x$cells)
		
		members <- membership[['Expressed']]
		non.members <- membership[['Remaining']]
		
		all.cells <- c(members, non.members)
		
		dist_matrix_cells <- dist_matrix[rownames(dist_matrix) %in% all.cells, colnames(dist_matrix) %in% all.cells]
		
		comparisons_members <- (length(members) * length(members) - 1) / 2
		comparisons_nonmembers <- (length(non.members) * length(non.members) - 1) /2
		comparisons_across <- length(members) * length(non.members)
		
		message(glue("number of comparisons are members = {comparisons_members}, non-members = {comparisons_nonmembers}, across = {comparisons_across}****\n"))
						
		groups <- ifelse(rownames(dist_matrix_cells) %in% members, "Members", "Non-members")
		distances_df <- dist_groups(dist_matrix_cells, groups)
		#distances_df <- distances_df[ ,c('Group1', 'Group2', 'Label', 'Distance')]
		distances_df <- distances_df[ ,c('Label', 'Distance')]
		distances_df$IC <- sig
		#distances_df
		distances_list[[sig]] <- distances_df	
		message(glue("done for {sig}***\n"))
    }
	#stopCluster(cl)
	#close(pb)
	message("combining dist data for all the components****")
	#distances_df <- do.call("rbind", distances_list)
	#distances_df <- do.call("rbind", distances_list)
	#distances_df$cellType <- cellType
	#distances_celltype_list[[cellType]] <- distances_df
	save(distances_list, file = glue("{outputDir}/distances_PC_singleCell_idhMut-{cellType}-{direction}Genes.Rda"))
	message(glue("done for {cellType}******\n\n"))
}





############obsolete#############
get.distances <- function(x, y, method = "euclidean") {
	combinations <- expand.grid(rownames(x), rownames(y))
	
	cl <- makeCluster(future::availableCores())
	registerDoParallel(cl)
	
	registerDoSNOW(cl)
	pb <- txtProgressBar(max = nrow(combinations), style = 3)
	progress <- function(n) setTxtProgressBar(pb, n)
	opts <- list(progress = progress)
	
	
    foreach(index = 1:nrow(combinations), .combine = "cbind", .options.snow = opts) %dopar% {
		cell1 <- combinations[index, 1]
		cell2 <- combinations[index, 2]
		
		data.dist <- rbind(x[cell1, ], y[cell2, ])
		dist(data.dist, method = method)
    }
	stopCluster(cl = NULL)
	close(pb)
	
	#distance.input <- rbind(x, y)
	#distances <- dist(distance.input, diag = FALSE, upper = FALSE)
	#return(distances)
}


#for (sig in names(meta.data)) {
#    membership <- meta.data %>% select(., select = paste(sig)) %>% rownames_to_column(., "cells")
#	membership <- split(membership, membership$select) %>% lapply(., function(x) x$cells)
#	
#	members <- membership[['Expressed']]
#	non.members <- membership[['Remaining']]

	#dist_matrix_members <- dist_subset(dist_matrix, members)
	#dist_matrix_across <- dist_subset(dist_matrix, members)
	
#	groups <- ifelse(rownames(dist_matrix) %in% members, "Members", "Non-members")
#	distances_df <- dist_groups(dist_matrix, groups)
#	distances_df <- distances_df[ ,c('Group1', 'Group2', 'Label', 'Distance')]
#	distances_df$IC <- sig
#	distances_list[[sig]] <- distances_df
	
	#dist_matrix_members <- dist_matrix[members, members]
	#dist_matrix_across <- dist_matrix[members, non.members]
	
	#distances_members <- dist_matrix_members[lower.tri(dist_matrix_members, diag = F)]
	#distances_across <- dist_matrix_across %>% data.frame %>% unlist()
	
	#distances_members_df <- data.frame(type = "members", distances = distances_members)
	#distances_across_df <- data.frame(type = "across", distances = distances_across)
   	#distances_df <- rbind(distances_members_df, distances_across_df)
	#distances_df$Dim <- sig
	#distances_df$direction <- direction
	#dim(distances_df)
	#distances_list[[sig]] <- distances_df
#	message(glue("done for {sig}***\n"))
#}