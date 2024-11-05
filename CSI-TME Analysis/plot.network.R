library(GGally)
library(network)
library(sna)
library(ggplot2)
library(igraph) ### to get adjacency matrix#####
library(glue)
library(dplyr)
library(RColorBrewer)
source('functions.network.R')


idh = "mut"
cohort <- "TCGA"
dims = "10"


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
#significant.interactions$prognosis <- ifelse(significant.interactions$HR > 0, "Worse", "Better")
#data.network <- significant.interactions[ ,c("partner1", "partner2", "bin", "prognosis")]

significant.interactions$type <- ifelse(significant.interactions$HR > 0, "Pro-tumor", "Anti-tumor")
data.network <- significant.interactions[ ,c("partner1", "partner2", "bin", "type")]


remove.bcells <- F

data.network$partner1 <- gsub("Dim.", "IC", data.network$partner1)
data.network$partner2 <- gsub("Dim.", "IC", data.network$partner2)



if (remove.bcells) {
	indexes <- grep("Bcell", significant.interactions$partners, invert = T)
	data.interactions <- significant.interactions[indexes, ]
	output.file <- glue("plots/cgin.components-IDH{idh}-tcga.cgga-technical-noBcells.svg")
} else {
	data.interactions <- significant.interactions
	output.file <- glue("plots/cgin.components-IDH{idh}-tcga.cgga-technical.svg")
}



#data.network <- data.interactions[ ,c("partner1", "partner2", "prognosis", "bin")] 

cgin.network <- create.network(data.network, multiple.arg = T, node.attr.names = c("cellType", "dimension"))
cgin.graph <- intergraph::asIgraph(cgin.network)

svglite::svglite(file = output.file, height = 12, width = 12)
plot.igraph.function(igraph = cgin.graph, col.nodes = "cellType", col.edges = "type", pallet.nodes = "Set2", pallet.edges = "Accent", position.legend.node = "topleft", position.legend.edge = "left", edge.size = 3, node.label.name = "dimension")
dev.off()  


cgin.degree <- degree(cgin.graph)
cellTypes <- V(cgin.graph)$cellType
dimension <- V(cgin.graph)$dimension

data.degree <- data.frame(cellType = cellTypes, dimension = dimension, degree = cgin.degree)

output.file <- glue("plots/degree.components.cellTypes.boxplot-IDH{idh}-tcga.cgga-technical.svg")
svg(file = output.file, height = 2, width = 4)
ggplot(data = data.degree, aes(x = cellType, y = degree)) + geom_boxplot() + theme_bw() + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 9))
dev.off()

#################interactions at the level of cell types######################
prognosis <- "Worse"
output.file <- glue("plots/cgin.cellTypes.{prognosis}Prognosis-IDH{idh}-tcga.cgga-technical.svg")

data.network <- significant.interactions[significant.interactions$prognosis == prognosis ,c("cellType1", "cellType2")] %>% magrittr::set_colnames(c("partner1", "partner2"))
data.network <- data.network %>% group_by(partner1, partner2) %>% tally() %>% data.frame()
cgin.network <- create.network(data.network, multiple.arg = T, node.attr.names = "cellType")
cgin.graph <- intergraph::asIgraph(cgin.network)
cgin.graph <- simplify(cgin.graph, edge.attr.comb = "sum")

svg(file = output.file, height = 4, width = 4)
plot.igraph.function(igraph = cgin.graph, col.nodes = "cellType", col.edges = "prognosis", pallet.nodes = "Set2", pallet.edges = "Accent", position.legend.node = "topleft", position.legend.edge = "left", edge.size = "n", node.label.name = "cellType")
dev.off()  

plot(cgin.graph, edge.width = E(cgin.graph)$n)

#####################plot LR network#################
idh = "mut"
cohort <- "TCGA"
dims = "10"


output.file <- glue("plots/cgin.LR-IDH{idh}-tcga.cgga-technical-noBcells.svg")

actual.file <- glue("ICA-JADE-{dims}.interactions.lr.filter-IDH{idh}-tcga.cgga-technical.txt")
lr.data <- read.table(actual.file, sep = "\t", header = T) %>% na.omit()
lr.data <- tidyr::separate(lr.data, col = communication, into = c("gene1", "gene2", "pairing"), sep = ":|-") 

partner1 <- paste(lr.data$cellType1, lr.data$gene1, sep = ":")
partner2 <- paste(l.data$cellType2, lr.data$gene2, sep = ":")
prognosis <- lr.data$prognosis
bin <- lr.data$bin
data.network <- data.frame(partner1, partner2, prognosis, bin)

lr.network <- create.network(data.network, multiple.arg = T, node.attr.names = c("cellType", "LR.gene"))
lr.graph <- intergraph::asIgraph(lr.network)

svg(file = output.file, height = 8, width = 8)
plot.igraph.function(igraph = lr.graph, col.nodes = "cellType", col.edges = "prognosis", node.label.name = "LR.gene", pallet.nodes = "Set2", pallet.edges = "Accent", position.legend.node = "topleft", position.legend.edge = "left")
dev.off()




#################LR at cell type level##################
idh = "mut"
cohort <- "TCGA"
dims = "10"

output.file <- glue("plots/cgin.LR-cellTypeLevel-IDH{idh}-tcga.cgga-technical.svg")

actual.file <- glue("ICA-JADE-{dims}.interactions.lr.filter-IDH{idh}-tcga.cgga-technical.txt")
lr.data <- read.table(actual.file, sep = "\t", header = T) %>% na.omit()

partner1 <- gsub(":\\w+", "", lr.data$lr.part1)
partner2 <- gsub(":\\w+", "", lr.data$lr.part2)
prognosis <- lr.data$prognosis
bin <- lr.data$bin
gene1 <- gsub("\\w+:", "", lr.data$lr.part1)
gene2 <- gsub("\\w+:", "", lr.data$lr.part2)
communication <- paste(gene1, gene2, sep = ":")
data.network <- data.frame(partner1, partner2, prognosis, bin, communication = communication)

lr.network <- create.network(data.network, multiple.arg = T, node.attr.names = c("cellType"))
lr.graph <- intergraph::asIgraph(lr.network)

svg(file = output.file, height = 8, width = 8)
plot.igraph.function(igraph = lr.graph, col.nodes = "cellType", col.edges = "prognosis", node.label.name = "cellType", edge.label.name = "communication", pallet.nodes = "Set2", pallet.edges = "Accent", position.legend.node = "topleft", position.legend.edge = "left")
dev.off()




################################################################
################################################################
################################################################
############plot subnetwork used for immunotherapy##############
################################################################
################################################################

actual.file <- glue("ICA-JADE-{dims}.FDR0.25.robustInteractions-IDH{idh}-tcga.cgga-technical.txt")

interactions <- read.table(actual.file, sep = "\t", header = T) %>% na.omit()
significant.interactions <- interactions %>% subset(., FDR.x < 0.20 & Reproduciblity > 7 & prognosis.x == prognosis.y & pvalue.y < 0.05)
significant.interactions <- significant.interactions[ ,c(2:7)]
colnames(significant.interactions) <- gsub("\\.x", "", colnames(significant.interactions))
significant.interactions <- format.interactions(significant.interactions)
significant.interactions <- get.partners(inputDf = significant.interactions)
#significant.interactions$prognosis <- ifelse(significant.interactions$HR > 0, "Worse", "Better")
#data.network <- significant.interactions[ ,c("partner1", "partner2", "bin", "prognosis")]

significant.interactions$type <- ifelse(significant.interactions$HR > 0, "Pro-tumor", "Anti-tumor")
data.network <- significant.interactions[ ,c("partner1", "partner2", "bin", "type")]


remove.bcells <- F

data.network$partner1 <- gsub("Dim.", "IC", data.network$partner1)
data.network$partner2 <- gsub("Dim.", "IC", data.network$partner2)



if (remove.bcells) {
	indexes <- grep("Bcell", significant.interactions$partners, invert = T)
	data.interactions <- significant.interactions[indexes, ]
	output.file <- glue("plots/cgin.components.subNetwork-IDH{idh}-tcga.cgga-technical-noBcells.svg")
} else {
	data.interactions <- significant.interactions
	output.file <- glue("plots/cgin.components.subNetwork-IDH{idh}-tcga.cgga-technical.svg")
}



#data.network <- data.interactions[ ,c("partner1", "partner2", "prognosis", "bin")] 

cgin.network <- create.network(data.network, multiple.arg = T, node.attr.names = c("cellType", "dimension"))
cgin.graph <- intergraph::asIgraph(cgin.network)

svglite::svglite(file = output.file, height = 8, width = 8)
plot.igraph.function(igraph = cgin.graph, col.nodes = "cellType", col.edges = "type", pallet.nodes = "Set2", pallet.edges = "Accent", position.legend.node = "topleft", position.legend.edge = "left", edge.size = 3, node.label.name = "dimension")
dev.off() 









##########################obsolete##################
#svg(file = output.file, height = 8, width = 8)
#ggnet2(cgin.network, color = "cellType", palette = "Set2", label = "dimension", label.size = 3.5)
#dev.off()

color.nodes <- brewer.pal(7,"Set2")
cellTypes <- V(cgin.graph)$cellType %>% unique()
names(color.nodes) <- cellTypes
V(cgin.graph)$color <- color.nodes[V(cgin.graph)$cellType]

color.edges <- brewer.pal(2,"Accent")[1:2]
prognosis <- E(cgin.graph)$prognosis %>% unique()
names(color.edges) <- prognosis
E(cgin.graph)$color <- color.edges[E(cgin.graph)$prognosis]

svg(file = output.file, height = 12, width = 12)
plot(cgin.graph, edge.curve = 0.1, 
	vertex.label = V(cgin.graph)$dimension, 
	vertex.color = V(cgin.graph)$color, 
	edge.color = E(cgin.graph)$color, 
	edge.width = 3, 
	vertex.size = 11, 
	layout = layout_with_fr)

	
legend("left", legend = names(color.nodes), fill = color.nodes, title = "celltype")
legend("bottomleft", legend = names(color.edges), fill = color.edges, title = "prognosis")
dev.off()  


########plot at the level of components########
plot.theme <- theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1), #face = "bold"),
	axis.text.y = element_text(size = 12, face = "bold"),
	strip.background =element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 12),
	#legend.title = element_blank(), legend.background = element_blank(), legend.position = "none",
	axis.title = element_text(size = 12), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	






















	##############Obsolete###############
cgin.data <- data.interactions[ ,c("partner1", "partner2", "prognosis")] 
cgin.adjacency <- graph_from_data_frame(cgin.data, directed = F) %>% as_adjacency_matrix(., type = "both") %>% as.matrix()
			
cgin.graph = network(cgin.adjacency, directed = FALSE, multiple = T)
cgin.vertices <- cgin.graph %v% "vertex.names"
cellTypes <- cgin.vertices %>% gsub(":.+", "", .)
dimensions <- cgin.graph %v% "vertex.names" %>% gsub(".+:", "", .)
dimensions <- gsub("Dim.", "IC", dimensions)

network::set.vertex.attribute(cgin.graph, attrname = "cellType", value = cellTypes)
network::set.vertex.attribute(cgin.graph, attrname = "dimension", value = dimensions)
network::set.edge.attribute(cgin.graph, attrname = "prognosis", value = dimensions)

svg(file = output.file, height = 8, width = 8)
ggnet2(cgin.graph, color = "cellType", palette = "Set2", label = dimensions, label.size = 3.5)
dev.off()

#set.network.attribute(cgin.graph, attrname = "cellType", value = cellTypes)
#set.edge.value(cgin.graph, attrname = "frequency", value = cgin.adjacency)
ggnet2(cgin.graph, color = "cellType")




############previous functions############

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

create.network <- function(data.network, direct.arg = F, multiple.arg = F, node.attr.names) {
	names.attributes <- colnames(data.network)[3:ncol(data.network)]
	count.attributes <- length(names.attributes)
	node.names <- c(data.network$partner1, data.network$partner2) %>% sort %>% unique()
	nodes <- 1:length(node.names)
	names(nodes) <- node.names
	
	data.network$head <- nodes[data.network$partner1]
	data.network$tail <- nodes[data.network$partner2]
	
	network.object <- network.initialize(length(nodes), directed = direct.arg, multiple = multiple.arg)
	network::add.edges(network.object, head = data.network$head, tail = data.network$tail)
	
	for (attribute in names.attributes) {
		network::set.edge.attribute(network.object, attrname = attribute, value = data.network[[attribute]])
	}

	#cellTypes <- gsub(":.+", "", node.names)
	#dimensions <- gsub(".+:", "", node.names)
	#dimensions <- gsub("Dim.", "IC", dimensions)
	
	node.attributes <- strsplit(node.names, ":")
	
	if (length(node.attr.names) > 0) {
		for (index in  1:length(node.attr.names)) {
			node.attr.name <- node.attr.names[index]
			node.attr.value <- lapply(node.attributes, function(x) x[index]) %>% unlist()
			network::set.vertex.attribute(network.object, attrname = node.attr.name, value = node.attr.value)
		}
	}
	
	
	#network::set.vertex.attribute(network.object, attrname = "cellType", value = cellTypes)
	#network::set.vertex.attribute(network.object, attrname = "dimension", value = dimensions)
	
	return(network.object)
}

plot.igraph.function <- function(igraph, col.nodes, col.edges, pallet.nodes, pallet.edges, node.label.name, edge.label.name, position.legend.node, position.legend.edge) {
	node.attr <- get.vertex.attribute(igraph, col.nodes) %>% unique()
	size.node.attr <- length(node.attr)
	
	color.nodes <- brewer.pal(size.node.attr, pallet.nodes)
	names(color.nodes) <- node.attr
	V(igraph)$color <- color.nodes[get.vertex.attribute(igraph, col.nodes)]
	
	edge.attr <- get.edge.attribute(igraph, col.edges) %>% unique()
	size.edge.attr <- length(edge.attr)
	color.edges <- brewer.pal(size.edge.attr,pallet.nodes)[1:size.edge.attr]
	names(color.edges) <- edge.attr
	E(igraph)$color <- color.edges[E(igraph)$prognosis]
	
	node.labels <- get.vertex.attribute(igraph, node.label.name)
	edge.labels <- get.edge.attribute(igraph, edge.label.name)
	
	plot(igraph, edge.curve = 0.1, 
			vertex.label = node.labels,
			edge.label = edge.labels,
			vertex.color = V(igraph)$color, 
			edge.color = E(igraph)$color, 
			edge.width = 3, 
			vertex.size = 11, 
			layout = layout_with_fr)
		
	legend(position.legend.node, legend = names(color.nodes), fill = color.nodes, title = col.nodes)
	legend(position.legend.edge, legend = names(color.edges), fill = color.edges, title = col.edges)
	return
}

plot.distribution <- function(igraph.object) {
	degree.values <- degree(igraph.object)
	degree.max <- max(degree.values)
	degree.distribution <- degree_distribution(igraph.object)
	degree.indexes <- c(0:degree.max)
	
	
	output <- fit_power_law(degree.distribution)
	
	output.plot <- plot(degree.indexes, degree.distribution)
	
	degree.data <- data.frame(degree.indexes, degree.distribution)
	
}

