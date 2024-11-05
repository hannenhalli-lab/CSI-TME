library(glue)
library(dplyr)




format.interactions <- function(interactions.df, input = "cellTypePairs") {
	if (input == "cellTypePairs") {
		tidyr::separate(interactions.df, col = cellTypes, into = c("cellType1", "cellType2"), sep = ":|-") %>% 
			tidyr::separate(., col = components, into = c("component1", "component2"), sep = ":|-") %>%
			mutate(., partner1 = paste(cellType1, component1, sep = ":"), partner2 = paste(cellType2, component2, sep = ":"))
	}
}

get.partners <- function(inputDf) {
	###must be outputed from format.interactions format####
	partners.df <- inputDf[ ,c('partner1', 'partner2')] %>% apply(., 1, sort) %>% t() %>% data.frame()
	colnames(partners.df) <- c('partner1', 'partner2')
	partners.df$partners <- paste(partners.df$partner1, partners.df$partner2, sep = "_")
	inputDf$partners <- partners.df$partners
	return(inputDf)
}

idh = "mut"
cohort <- "TCGA"
dims = "10"
database.source <- "cellchatdb"

############################################################################
############################################################################
############get adjacency matrix from ligand receptor data############

actual.file <- glue("ICA-JADE-{dims}.FDR0.25.robustInteractions-IDH{idh}-tcga.cgga-technical.txt")

interactions <- read.table(actual.file, sep = "\t", header = T) %>% na.omit()
significant.interactions <- interactions %>% subset(., FDR.x < 0.20 & Reproduciblity > 7)
significant.interactions <- significant.interactions[ ,c(2:7)]
colnames(significant.interactions) <- gsub("\\.x", "", colnames(significant.interactions))
significant.interactions <- format.interactions(significant.interactions)
significant.interactions <- get.partners(inputDf = significant.interactions)
significant.interactions$prognosis <- ifelse(significant.interactions$HR > 0, "Worse", "Better")

lr.files <- list.files("./", pattern = glue("ligand.receptor.{database.source}")) %>% grep("\\.txt", ., value = T)

lr.files.data <- lapply(lr.files, read.table, sep = "\t", header = T)
lr.df <- bind_rows(lr.files.data)
#lr.df <- read.table(glue("ligand.receptor.{database.source}-ICA-JADE-{dims}.signatureGenes.{cohort}.positiveGenes.txt"), sep = "\t", header = T)

lr.df <- format.interactions(lr.df)
lr.df <- get.partners(inputDf = lr.df) #####the interacting pairs are sorted alphbetically and loose their ordering in the variable named as 'partners'
lr.df$communication <- paste(lr.df$Gene1, lr.df$Gene2, lr.df$type, sep = ":")
lr.df$lr.part1 <- paste(lr.df$cellType1, lr.df$Gene1, sep = ":")
lr.df$lr.part2 <- paste(lr.df$cellType2, lr.df$Gene2, sep = ":")

lr.df <- lr.df[ ,c("partners", "lr.part1", "lr.part2", "communication", 'direction.genes')]


data.network.lr <- merge(significant.interactions, lr.df, by = "partners")
write.table(data.network.lr, file = glue("ICA-JADE-{dims}.interactions.lr.filter-IDH{idh}-tcga.cgga-technical.txt"), sep = "\t", quote = F, row.names = F)


distribution.interactions <- table(significant.interactions$bin, significant.interactions$prognosis)
distribution.lr <- table(data.network.lr$bin, data.network.lr$prognosis)

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################


###########################for random interactions#########################

lr.df <- read.table(glue("ligand.receptor.{database.source}-ICA-JADE-{dims}.signatureGenes.{cohort}.positiveGenes.txt"), sep = "\t", header = T)
lr.df <- format.interactions(lr.df)
lr.df <- get.partners(inputDf = lr.df)
lr.df$communication <- paste(lr.df$Gene1, lr.df$Gene2, lr.df$type, sep = ":")
lr.df <- lr.df[ ,c("partners", "communication")]



random.file <- "ICA-JADE-10.randomInteractions-IDHmut-TCGA-technical.txt"
interactions <- read.table(random.file, sep = "\t", header = T) %>% na.omit()
significant.interactions <- interactions %>% subset(., FDR < 0.25)
significant.interactions <- format.interactions(significant.interactions)
significant.interactions <- get.partners(inputDf = significant.interactions)
significant.interactions$prognosis <- ifelse(significant.interactions$HR > 0, "Worse", "Better")


data.random.lr.list <- list()
for (iteration in 1:10) {
	significant.iteration <- significant.interactions[significant.interactions$iteration == iteration, ]
	data.random.lr <- merge(significant.iteration, lr.df, by = "partners")
	data.random.lr.list[[iteration]] <- data.random.lr$partners %>% unique()
}
















