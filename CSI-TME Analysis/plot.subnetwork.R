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


actual.file <- glue("ICA-JADE-{dims}.FDR0.25.robustInteractions-IDH{idh}-tcga.cgga-technical.txt")

interactions <- read.table(actual.file, sep = "\t", header = T) %>% na.omit()
significant.interactions <- interactions %>% subset(., FDR.x < 0.20 & Reproduciblity > 7)
significant.interactions <- significant.interactions[ ,c(1:7)]
colnames(significant.interactions) <- gsub("\\.x", "", colnames(significant.interactions))
significant.interactions <- format.interactions(significant.interactions)
significant.interactions <- get.partners(inputDf = significant.interactions)
significant.interactions$prognosis <- ifelse(significant.interactions$HR > 0, "Pro-Tumor", "Anti-Tumor")


#filter.interaction <- c("Malignant:Dim.5", "Malignant:Dim.6", "Malignant:Dim.7")
filter.interaction <- c("Malignant:Dim.5", "Malignant:Dim.6", "Malignant:Dim.7")
#output.file <- glue("plots/cgin.{filter.interaction}-IDH{idh}-tcga.cgga-technical.svg")
output.file <- glue("plots/cgin.Stemness-IDH{idh}-tcga.cgga-technical.svg")


data.network <- significant.interactions[ ,c("partner1", "partner2", "bin", "prognosis")] %>% .[.$partner1 %in% filter.interaction | .$partner2 %in% filter.interaction, ]

data.network$partner1 <- gsub("Dim.", "IC", data.network$partner1)
data.network$partner2 <- gsub("Dim.", "IC", data.network$partner2)


cgin.network <- create.network(data.network, multiple.arg = T, direct.arg = T, node.attr.names = c("cellType", "dimension"))
cgin.graph <- intergraph::asIgraph(cgin.network)



svglite::svglite(file = output.file, height = 6, width = 6)
plot.igraph.function(igraph = cgin.graph, col.nodes = "cellType", col.edges = "prognosis", 
		node.label.name = "dimension", edge.label.name = "bin", 
		pallet.nodes = "Set2", pallet.edges = "Accent",
		arrow.shape.name = "bin", arrow.shape.criteria = "Bin.3", 
		position.legend.node = "topleft", position.legend.edge = "left")
dev.off()  




###################based on mutations#####################
mutation.interactions <- read.table(glue("ICA-JADE-{dims}.FDR0.25.robustInteractions-IDH{idh}-association.Mutations.txt"), sep = "\t", header = T)

mutation.interactions$Association <- ifelse(mutation.interactions$Estimate > 0 , "Positive", "Negative")

association.significant <- mutation.interactions[mutation.interactions$FDR < 0.25 & mutation.interactions$Estimate > 1,  ]
association.significant$significance <- -log10(association.significant$FDR)

association.significant <- association.significant[ ,c('predictor', 'response')]
association.significant <- association.significant[!(association.significant$response %in% c("TP53", "ATRX")), ]
association.significant <- association.significant %>% group_by(predictor) %>% summarise(mutations = paste(response,collapse=','))

data.network <- merge(significant.interactions, association.significant, by.x = "string", by.y = "predictor")
data.network <- data.network[ ,c("partner1", "partner2", "bin", "prognosis", "mutations")]

data.network$partner1 <- gsub("Dim.", "IC", data.network$partner1)
data.network$partner2 <- gsub("Dim.", "IC", data.network$partner2)

cgin.network <- create.network(data.network, multiple.arg = T, direct.arg = T, node.attr.names = c("type", "dimension"))
cgin.graph <- intergraph::asIgraph(cgin.network)


output.file <- glue("plots/cgin.mutations-IDH{idh}-tcga.cgga-technical.svg")

svg(file = output.file, height = 9, width = 9)
plot.igraph.function(igraph = cgin.graph, col.nodes = "cellType", col.edges = "prognosis", 
		node.label.name = "dimension", edge.label.name = "mutations", 
		pallet.nodes = "Set2", pallet.edges = "Accent", node.size = 8,
		arrow.shape.name = NULL, arrow.shape.criteria = "Bin.3", layout.name = layout_as_bipartite,
		position.legend.node = "topleft", position.legend.edge = "left")
dev.off()  








###############plot for genes###########

cellType1 <- "Malignant"
cellType2 <- "Tcell"

component1 <- "Dim.7"
component2 <- "Dim.7"

direction1 <- "negative"
direction2 <- "negative"

covariateType <- "technical"

name.list <- glue("{cellType1}.{component1}.{direction1}-{cellType2}.{component2}.{direction2}")
input.file <- glue("{name.list}-IDH{idh}-{cohort}-{covariateType}.txt")

interactions <- read.table(input.file, sep = "\t", header = T) %>% na.omit()
interactions <- interactions[interactions$bin == "Bin.1" & interactions$FDR < 0.20,]
interactions$prognosis <- ifelse(interactions$HR > 0, "Pro-tumor", "Anti-tumor")

data.network <- interactions[ ,c("cellTypes", "components", "prognosis")]
data.network <- format.interactions(data.network)

cgin.network <- create.network(data.network, multiple.arg = T, direct.arg = T, node.attr.names = c("cellType", "dimension"))
cgin.graph <- intergraph::asIgraph(cgin.network)

output.file <- glue("plots/network.{name.list}-IDH{idh}-{cohort}-{covariateType}.svg")
svg(file = output.file, height = 14, width = 14)
plot.igraph.function(igraph = cgin.graph, col.nodes = "type", col.edges = "prognosis", 
		node.label.name = "dimension", edge.label.name = "bin", 
		pallet.nodes = "Set2", pallet.edges = "Accent",
		arrow.shape.name = NULL, arrow.shape.criteria = NA, 
		position.legend.node = "topleft", position.legend.edge = "left")
dev.off()  


#########################plot ligand receptor network######################

covariateType <- "technical"; idh = "mut"; cohort <- "tcga"
lr.file <- "ICA-JADE-10.interactions.lr.filter.simultaneous.Status-IDHmut-tcga.cgga-technical.txt"

interactions <- read.table(lr.file, sep = "\t", header = T)
interactions$type <- ifelse(interactions$prognosis == "Better", "Anti-tumor", "Pro-tumor")

#data.network <- interactions[ ,c("lr.part1", "lr.part2", "bin", "type")]
#names(data.network) <- c("partner1", "partner2", "bin", "type")

data.network <- interactions[ ,c("lr.part1", "lr.part2", "lr.activity", "type")]
names(data.network) <- c("partner1", "partner2", "lr.activity", "type")


cgin.network <- create.network(data.network, multiple.arg = T, direct.arg = T, node.attr.names = c("cellType", "Gene"))
cgin.graph <- intergraph::asIgraph(cgin.network)
E(cgin.graph)$lty <- ifelse(E(cgin.graph)$lr.activity == "Active", 1, 3)

name.list <- "lr.interactions"
output.file <- glue("plots/network.{name.list}-IDH{idh}-{cohort}-{covariateType}.svg")
svglite::svglite(file = output.file, height = 15, width = 15)
par(mar=c(0,0,0,0))
plot.igraph.function(igraph = cgin.graph, col.nodes = "cellType", col.edges = "type", 
		node.label.name = "Gene", node.label.size = 0.8, node.size = 8,
		edge.label.name = "bin", edge.label.size = 0.8,
		pallet.nodes = "Set2", pallet.edges = "Accent",
		arrow.shape.name = NULL, arrow.shape.criteria = NA, 
		position.legend.node = "topleft", position.legend.edge = "left")
dev.off()  


function(igraph, col.nodes, col.edges, pallet.nodes, pallet.edges, node.label.name, edge.label.name = NULL, arrow.shape.name = NULL, arrow.shape.criteria = NA, position.legend.node, position.legend.edge, edge.size = 3) {
	##igraph --> an igraph object
	##col.nodes --> node attribute used to color nodes
	##col.edges --> edge attribute used to color edges
	##pallet.nodes --> color pallet for nodes
	##pallet.edges --> color pallet for edges
	##node.label.name --> node attribute used to print at nodes in graph
	##edge.label.name --> edge attribute used to print at edges in graph
	##arrow.shape.name --> edge attribute used to shape the edges in graph
	##arrow.shape.criteria --> one of the value from arrow.shape.name in graph to select shapes
	##position.legend.node --> where to position the node legend
	##position.legend.edge --> where to position the edge legend
	##edge.size --> either a numbner corresponding to width of edge or an edge attribute based on which to size the edges
	
	node.attr <- get.vertex.attribute(igraph, col.nodes) %>% unique()
	size.node.attr <- length(node.attr)
	
	color.nodes <- brewer.pal(size.node.attr, pallet.nodes)
	color.nodes <- color.nodes[1:size.node.attr]
	names(color.nodes) <- node.attr
	V(igraph)$color <- color.nodes[get.vertex.attribute(igraph, col.nodes)]
	
	edge.attr <- get.edge.attribute(igraph, col.edges) %>% unique()
	size.edge.attr <- length(edge.attr)
	color.edges <- brewer.pal(size.edge.attr,pallet.nodes)[1:size.edge.attr]
	names(color.edges) <- edge.attr
	E(igraph)$color <- color.edges[get.edge.attribute(igraph, col.edges)]
	
	node.labels <- get.vertex.attribute(igraph, node.label.name)
	
	if (is.numeric(edge.size)) {
		edge.size <- edge.size
	}
	if (is.character(edge.size)) {
		edge.size <- get.edge.attribute(igraph, edge.size)
	}
	
	if (!is.null(arrow.shape.name)) {
		arrow.values <- get.edge.attribute(igraph, arrow.shape.name)
		arrow.shape <- ifelse(arrow.values %in% arrow.shape.criteria, ">", "-")
	} else {
		arrow.shape <- rep("-", E(igraph) %>% length)
	}
	
	if (layout.name == "layout_as_bipartite") {
		types.arg <- node.attr[1] == get.vertex.attribute(igraph, "type")
		igraph = layout_as_bipartite(graph = igraph, types = types.arg)
	}
	if (!is.null(edge.label.name)) {
		edge.labels <- get.edge.attribute(igraph, edge.label.name)
		plot(igraph, edge.curve = 0.1, 
				vertex.label = node.labels,
				vertex.label.cex = node.label.size,
				edge.label = edge.labels,
				vertex.color = V(igraph)$color, 
				edge.color = E(igraph)$color,
				edge.width = edge.size, 
				vertex.size = 11,
				edge.arrow.mode = arrow.shape,
				types = types.arg,
				layout = get(layout.name))
	} else {
		plot(igraph, edge.curve = 0.1, 
				vertex.label = node.labels,
				vertex.label.cex = node.label.size,
				vertex.color = V(igraph)$color, 
				edge.color = E(igraph)$color, 
				edge.width = edge.size, 
				vertex.size = 11, 
				layout = get(layout.name))
	}
	
		
	legend(position.legend.node, legend = names(color.nodes), fill = color.nodes, title = col.nodes)
	legend(position.legend.edge, legend = names(color.edges), fill = color.edges, title = col.edges)
	return