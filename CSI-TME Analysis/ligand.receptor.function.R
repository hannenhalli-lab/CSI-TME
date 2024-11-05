library(dplyr); library(magrittr); library(glue)
library(doParallel); library(tidyr)
source('../functionsForInteractionAnalysis.R')

format.interactions <- function(interactions.df, input = "cellTypePairs") {
	if (input == "cellTypePairs") {
		tidyr::separate(interactions.df, col = cellTypes, into = c("cellType1", "cellType2"), sep = ":|-") %>% 
			tidyr::separate(., col = components, into = c("component1", "component2"), sep = ":|-") %>%
			mutate(., partner1 = paste(cellType1, component1, sep = ":"), partner2 = paste(cellType2, component2, sep = ":"))
	}
}

database.source <- "cellchatdb"
if (database.source == "cellchatdb") {
	load('~/Downloads/CellChatDB.human.rda')
	database <- CellChatDB.human[['interaction']] %>% subset(., select = c(ligand, receptor)) %>% unique() %>% set_colnames(c("Ligand", "Receptor"))
	database <- tidyr::separate_rows(database, Ligand, Receptor, sep = "_", convert = FALSE) %>% mutate_all(., .funs = toupper) %>% unique()
	
}



idh = "mut"
cohort <- "TCGA"
dims = "10"
load(glue("ICA-JADE-{dims}.signatureGenes.{tolower(cohort)}.IDH{idh}--filtered-noNorm-optimalRank.Rda"))

cellTypes <- names(signature.list)
cellTypeCombinations <- combn(cellTypes, 2)

ligand.receptor.interaction.list <- list(); iteraction.iterations <- list(); 
#direction.genes = "positive:positive"; 
#direction.genes = "negative:positive"; 
#direction.genes = "positive:negative"; 
direction.genes = "negative:negative"; 

#direction.genes <- "np"
NM_list <- list()
NM_size <- 0
for (cellTypeIndex in 1:ncol(cellTypeCombinations)) {
	
		cellType1 <- cellTypeCombinations[1, cellTypeIndex]
		cellType2 <- cellTypeCombinations[2, cellTypeIndex]
		
		cellType1.sig <- signature.list[[cellType1]]
		cellType2.sig <- signature.list[[cellType2]]
		
		lv.1 <- names(cellType1.sig); lv.2 <- names(cellType2.sig)
		componentCombinations <- expand.grid(lv.1, lv.2)
		
		sizeCount <- nrow(componentCombinations)
		cores=detectCores()
		cl <- makeCluster(cores[1] - 2) #not to overload your computer
		registerDoParallel(cl)
		
		componentInteractionList <- NULL
		#NM_size <- 0
		componentInteractionList <- foreach(componentIndex=1:sizeCount, .packages=c('dplyr', 'magrittr')) %dopar% 
			{
				component1 <- componentCombinations[componentIndex, 1]
				component2 <- componentCombinations[componentIndex, 2]
				
				directions <- strsplit(direction.genes, ":") %>% unlist()
				
				direction1.genes <- directions[1];
				direction2.genes <- directions[2];
				
				
				component1.genes <- cellType1.sig[[component1]][[direction1.genes]]
				component2.genes <- cellType2.sig[[component2]][[direction2.genes]]
				
				
				message(NM_size)
				
				if (length(component1.genes) > 2 & length(component2.genes) > 2) {
					#lr.df <- get.LR.interactions(list1 = component1.genes, list2 = component2.genes, database = database)
					#lr.df
					NM <- length(component1.genes) * length(component2.genes)
					#NM_size <- NM_size + NM
					NM
				} else {
					NA
				}

			}
			stopCluster(cl)
			componentNames <- glue("{componentCombinations[, 1]}:{componentCombinations[, 2]}")
			names(componentInteractionList) <- componentNames
			ligand.receptor.interaction.list[[cellTypeIndex]] <- componentInteractionList
			print(glue("done for {cellTypeIndex} cellType combinations"))
}
names(ligand.receptor.interaction.list) <- glue("{cellTypeCombinations[1, ]}:{cellTypeCombinations[2, ]}")
save(ligand.receptor.interaction.list, file = glue("ligand.receptor.{database.source}-ICA-JADE-{dims}.signatureGenes.{cohort}.{direction.genes}.Genes.Rda"))

lr.df <- lapply(ligand.receptor.interaction.list, function(x) {
	x <- x[!is.na(x)];
	bind_rows(x, .id = "components")
}) %>% bind_rows(., .id = "cellTypes")

lr.df <- na.omit(lr.df)

lr.df$direction.genes <- direction.genes
write.table(lr.df, file = glue("ligand.receptor.{database.source}-ICA-JADE-{dims}.signatureGenes.{cohort}.{direction.genes}.Genes.txt"), sep = "\t", quote = F, row.names = F)


#####get sum of N*M values###
NM <- unlist(unlist(ligand.receptor.interaction.list))





#######################mostly obsolete#######################



############################################################################
############################################################################
############get adjacency matrix from ligand receptor data############
library(igraph)
lr.df <- read.table(glue("ligand.receptor.{database.source}-ICA-JADE-{dims}.signatureGenes.{cohort}.positiveGenes.txt"), sep = "\t", header = T)
lr.df <- format.interactions(lr.df)

#lr.df %>% group_by(cellTypes, components, type, size.ligands, size.receptors) %>% tally
lr.degree <- lr.df %>% group_by(partner1, partner2) %>% tally()

lr.adjacency <- lr.degree %>% graph_from_data_frame(., directed = F) %>% 
	 	as_adjacency_matrix(., type = "both", attr = "n") %>% as.matrix()


actual.file <- glue("ICA-JADE-{dims}.FDR0.25.robustInteractions-IDH{idh}-tcga.cgga-technical.txt")
random.file <- glue("ICA-JADE-{dims}.randomInteractions-IDH{idh}-TCGA-technical.txt")

interactions <- read.table(actual.file, sep = "\t", header = T) %>% na.omit()
significant.interactions <- interactions %>% subset(., FDR.x < 0.25 & Reproduciblity > 7 & pvalue.y < 0.2 & prognosis.x == prognosis.y )
significant.interactions <- significant.interactions[ ,c(2:7)]
colnames(significant.interactions) <- gsub("\\.x", "", colnames(significant.interactions))

random <- read.table(random.file, sep = "\t", header = T) %>% na.omit()
random.interactions <- random[random$FDR < 0.20, ]

significant.interactions <- format.interactions(significant.interactions)
random.interactions <- format.interactions(random.interactions)



lr.actual <- list()
for (row.index in 1:nrow(significant.interactions)) {
	partner1 <- significant.interactions[row.index, 'partner1']
	partner2 <- significant.interactions[row.index, 'partner2']
	
	lr <- try(lr.adjacency[partner1, partner2], silent = T)
	if (!is.numeric(lr)) {
		lr <- 0
	}
	lr.actual[[row.index]] <- lr
}
lr.random <- list()

for (row.index in 1:nrow(random.interactions)) {
	partner1 <- random.interactions[row.index, 'partner1']
	partner2 <- random.interactions[row.index, 'partner2']
	
	lr <- try(lr.adjacency[partner1, partner2], silent = T)
	if (!is.numeric(lr)) {
		lr <- 0
	}
	lr.random[[row.index]] <- lr
}

significant.interactions$lr <- unlist(lr.actual)
random.interactions$lr <- unlist(lr.random)
boxplot(random.interactions$lr, significant.interactions$lr)

##################################################################################################################
##################################################################################################################
##################################################################################################################
lr.df <- read.table(glue("ligand.receptor.{database.source}-ICA-JADE-{dims}.signatureGenes.{cohort}.positiveGenes.txt"), sep = "\t", header = T)

lr.df <- lr.df %>% group_by(cellTypes, components, type, size.ligands, size.receptors) %>% tally
lr.df <- format.interactions(lr.df)
lr.df$string <- paste(lr.df$partner1, lr.df$partner2, sep = "_")

lr.total <- nrow(database)
database.list <- split(database, database$Ligand)

ligands.all <- names(database.list) %>% unique
receptor.all <- lapply(database.list, function(x) x$Receptor) %>% unlist() %>% unique()

ligands.total <- length(ligands.all)
receptor.total <- length(receptor.all)
global.expectation <- lr.total / (ligands.total * receptor.total)


####Method 1###### 
#observed <- lr.df$n / (lr.df$size.ligands * lr.df$size.receptors)
#expected <- lr.total / (ligands.total * receptor.total)
#lr.df$OE <- observed / expected

####Method 2######both methods are same####
observed <- lr.df$n
expected <- lr.df$size.ligands * lr.df$size.receptors * global.expectation
lr.df$OE <- observed / expected


#interaction.file <- glue("ICA-JADE-{dims}.FDR0.25.robustInteractions-IDH{idh}-tcga.cgga-technical.txt")
interaction.file <- "ICA-JADE-10.allInteractions-IDHmut-TCGA-technical.txt"
interaction.data <- read.table(interaction.file, sep = "\t", header = T)
interaction.data$group <- ifelse(interaction.data$FDR < 0.25, "interaction", "background")

interaction.data <- format.interactions(interaction.data)

interaction.sorted <- interaction.data[ ,c('partner1', 'partner2', 'group')]
interaction.sorted <- interaction.sorted[ ,c('partner1', 'partner2')] %>% apply(., 1, sort) %>% t()
interaction.sorted <- data.frame(interaction.sorted) %>% set_colnames(c('partner1', 'partner2'))
interaction.sorted$group <- interaction.data$group
interaction.sorted$string <- paste(interaction.sorted$partner1, interaction.sorted$partner2, sep = "_")
interaction.sorted$bin <- interaction.data$bin

interactions.lr <- merge(interaction.sorted, lr.df, by = "string", all.x = T)
interactions.lr$OE <- ifelse(is.na(interactions.lr), 0, interactions.lr$OE)
#boxplot(OE ~ group, data = interactions.lr, outline = F)
ggplot(data = interactions.lr, aes(x = group, y = OE)) + geom_boxplot() + facet_wrap(~bin, scale = "free_y")




