library(dplyr); library(magrittr); library(glue)
library(tidyr); library(igraph)
library(ggplot2); library(ggpubr)
source('../functionsForInteractionAnalysis.R')

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

database.source <- "cellchatdb"
if (database.source == "cellchatdb") {
	load('~/Downloads/CellChatDB.human.rda')
	database <- CellChatDB.human[['interaction']] %>% subset(., select = c(ligand, receptor)) %>% unique() %>% set_colnames(c("Ligand", "Receptor"))
	database <- tidyr::separate_rows(database, Ligand, Receptor, sep = "_", convert = FALSE) %>% mutate_all(., .funs = toupper) %>% unique()
}


idh = "mut"
cohort <- "TCGA"
dims = "10"
database.source <- "cellchatdb"
onlyInteractions <- T

remove.bcells <- F

########load interactions and get degree of interaction between cell types#########
#interaction.file <- glue("ICA-JADE-{dims}.allInteractions-IDH{idh}-TCGA-technical.txt")
interaction.file <- glue("ICA-JADE-{dims}.FDR0.25.robustInteractions-IDH{idh}-tcga.cgga-technical.txt")
interaction.data <- read.table(interaction.file, sep = "\t", header = T)
interaction.data <- interaction.data %>% subset(., FDR.x < 0.20 & Reproduciblity > 7)
#interaction.data <- interaction.data %>% subset(., pvalue < 0.05)
#interaction.data <- interaction.data %>% subset(., FDR.x < 0.20 & Reproduciblity > 7 & pvalue.y < 0.2 & prognosis.x == prognosis.y )


if (remove.bcells) {
	indexes <- grep("Bcell", interaction.data$string, invert = T)
	interaction.data <- interaction.data[indexes, ]
}

interaction.data <- interaction.data[ ,c('cellTypes.x', 'components.x')] %>% set_colnames(c("cellTypes", "components"))

interaction.data <- format.interactions(interaction.data)
interaction.data <- get.partners(inputDf = interaction.data)



#######load LR pairs among the ICs of cell type and get mean OE for cell type pairs######
#direction.genes <- "positive"

direction.genes <- c("positive:positive", "negative:negative", "negative:positive", "positive:negative")
lr.df.list <- list()
for (direction in direction.genes) {
	lr.df <- read.table(glue("ligand.receptor.{database.source}-ICA-JADE-{dims}.signatureGenes.{cohort}.{direction}.Genes.txt"), sep = "\t", header = T)
	message(nrow(lr.df))
	#######count LR and RL pairs for groups
	lr.df <- lr.df %>% group_by(cellTypes, components, type, size.ligands, size.receptors) %>% tally()
	lr.df$string <- paste(lr.df$cellTypes, lr.df$components, sep = "_")
	lr.df.list[[direction]] <- lr.df
}
lr.df <- bind_rows(lr.df.list, .id = "direction")
lr.df.combined <- lr.df %>% group_by(cellTypes, components, type, string) %>% summarize(size.ligands = sum(size.ligands), size.receptors = sum(size.receptors), n = sum(n))
lr.df.combined$direction = "combined"

lr.df.combined <- lr.df.combined[ ,sort(colnames(lr.df.combined))]
lr.df <- lr.df[ ,sort(colnames(lr.df))]


lr.df <- rbind(lr.df, lr.df.combined)


if (remove.bcells) {
	indexes <- grep("Bcell", lr.df$string, invert = T)
	lr.df <- lr.df[indexes, ]
}


lr.total <- nrow(database)
database.list <- split(database, database$Ligand)

ligands.all <- names(database.list) %>% unique
receptor.all <- lapply(database.list, function(x) x$Receptor) %>% unlist() %>% unique()

ligands.total <- length(ligands.all)
receptor.total <- length(receptor.all)
global.expectation <- lr.total / (ligands.total * receptor.total)



observed <- lr.df$n
expected <- lr.df$size.ligands * lr.df$size.receptors * global.expectation
lr.df$OE <- observed / expected

lr.df <- format.interactions(lr.df)
lr.df <- get.partners(lr.df)


#####following two lines to only retain the LR pairs among the interacting components#####
if(onlyInteractions == T) {
	lr.df <- lr.df[lr.df$partners %in% interaction.data$partners, ] #### select LRs among only the interacting components
	interaction.sorted <- interaction.data[interaction.data$partners %in%  lr.df$partners, c('cellType1', 'cellType2')] %>% apply(., 1, sort) %>% t()
	
}  else {
	interaction.sorted <- interaction.data[ ,c('cellType1', 'cellType2')] %>% apply(., 1, sort) %>% t()
}

############
interaction.sorted <- data.frame(interaction.sorted) %>% set_colnames(c('cellType1', 'cellType2'))
interaction.degree <- interaction.sorted %>% group_by(cellType1, cellType2) %>% tally()
interaction.degree$pair <- paste(interaction.degree$cellType1, interaction.degree$cellType2, sep = ":")




lr.meanOE <- lr.df %>% group_by(cellType1, cellType2, type, direction) %>% summarize(OE = mean(OE))
lr.meanOE <- lr.meanOE %>% group_by(cellType1, cellType2, direction) %>% summarize(OE = mean(OE))

#lr.meanOE <- lr.df %>% group_by(cellType1, cellType2, type) %>% summarize(OE = sum(n))
#lr.meanOE <- lr.df %>% group_by(cellType1, cellType2) %>% summarize(OE = sum(OE))
lr.meanOE$pair <- paste(lr.meanOE$cellType1, lr.meanOE$cellType2, sep = ":")

ggTheme <- theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5), #face = "bold"),
	axis.text.y = element_text(size = 12, face = "bold"),
	strip.background =element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 12),
	#legend.title = element_blank(), legend.background = element_blank(), legend.position = "none",
	axis.title = element_text(size = 12), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	
	

data.plot <- merge(interaction.degree, lr.meanOE, by = "pair")
scatter.plot <- ggplot(data = data.plot, aes(x = n, y =  OE))  + theme_bw() + facet_wrap(~direction, scale = "free", nrow = 1) +
				geom_point(aes(size = 3)) + geom_smooth(method = "lm", se = F) + ggpubr::stat_cor(method = "pearson", col = "red") +
				xlab("genetic interaction degree") + ylab("Observed by Expected of Ligand Receptors") 

svglite::svglite(glue("plots/ligand.receptor.scatterPlot-allSignaturePairs.svg"), height = 3, width = 10)
scatter.plot + ggTheme
dev.off()


####for main figure###
data.plot <- data.plot[data.plot$direction == "combined", ]

scatter.plot <- ggplot(data = data.plot, aes(x = n, y =  OE))  + theme_bw() + facet_wrap(~direction, scale = "free", nrow = 1) +
				geom_point(aes(size = 3)) + geom_smooth(method = "lm", se = F) + ggpubr::stat_cor(method = "pearson", col = "red") +
				xlab("genetic interaction degree") + ylab("Observed by Expected of Ligand Receptors") 

svglite::svglite(glue("plots/ligand.receptor.scatterPlot-combined.svg"), height = 3.5, width = 3.5)
scatter.plot + ggTheme
dev.off()

if (onlyInteractions == T) {
	svglite::svglite(glue("plots/ligand.receptor.scatterPlot-OnlyInteractions.svg"), height = 4, width = 4)
	scatter.plot + ggTheme
	dev.off()
	
	
	data.plot$level <- factor(data.plot$n)
	box.plot <- ggplot(data = data.plot, aes(x = level, y =  OE))  + theme_bw() + 
					geom_boxplot() + ggpubr::stat_compare_means(label = "p.format") +
					xlab("degree LR") + ylab("Observed by Expected of Ligand Receptors") 
					
	svglite::svglite(glue("plots/ligand.receptor.boxplot-OnlyInteractions.svg"), height = 3.5, width = 2) 
	box.plot + ggTheme
	dev.off()
}
if (onlyInteractions == F) {
	svglite::svglite(glue("plots/ligand.receptor.scatterPlot.svg"), height = 4, width = 4)
	scatter.plot + ggTheme
	dev.off()
}











	###############Using randomized interaction data##################
random.file <- glue("ICA-JADE-{dims}.randomInteractions-IDH{idh}-TCGA-technical.txt")
random.data <- read.table(random.file, sep = "\t", header = T) %>% na.omit()
random.data <- random.data[random.data$FDR < 0.30, ]

random.data <- format.interactions(random.data)

random.sorted <- random.data[ ,c('cellType1', 'cellType2')] %>% apply(., 1, sort) %>% t()
random.sorted <- data.frame(random.sorted) %>% set_colnames(c('cellType1', 'cellType2'))
random.sorted$iteration <- random.data$iteration

random.degree <- random.sorted %>% group_by(cellType1, cellType2, iteration) %>% tally()
random.degree$pair <- paste(random.degree$cellType1, random.degree$cellType2, sep = ":")


data.plot <- merge(random.degree, lr.meanOE, by = "pair")
ggplot(data = data.plot, aes(x = n, y =  OE))  + theme_bw() + facet_wrap(~iteration, scale = "free") +
	geom_point(aes(size = 3)) + geom_smooth(method = "lm", se = F) + ggpubr::stat_cor(method = "spearman", col = "red") +
	xlab("genetic interaction degree") + ylab("Observed by Expected of Ligand Receptors") 




