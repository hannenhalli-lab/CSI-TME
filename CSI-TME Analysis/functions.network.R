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

create.network <- function(data.network, direct.arg = F, multiple.arg = F, node.attr.names = NULL) {
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

	cellTypes <- gsub(":.+", "", node.names)
	dimensions <- gsub(".+:", "", node.names)
	dimensions <- gsub("Dim.", "IC", dimensions)
	
	node.attributes <- strsplit(node.names, ":")
	
	if (!is.null(node.attr.names)) {
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

plot.igraph.function <- function(igraph, col.nodes, col.edges, pallet.nodes, pallet.edges, node.label.name, edge.label.name = NULL, arrow.shape.name = NULL, arrow.shape.criteria = NA, position.legend.node, position.legend.edge, edge.size = 3, node.size = 8, node.label.size = 1, edge.label.size = 1) {
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
	##node.size --> a number indicating size of the node plotted
	##node.label.size --> a number indicating the size of fonts for node label
	##edge.label.size --> a number indicating the size of fonts for edge label
	
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
	
	if (!is.null(edge.label.name)) {
		edge.labels <- get.edge.attribute(igraph, edge.label.name)
		plot(igraph, edge.curve = 0.1, 
				vertex.label = node.labels,
				vertex.label.cex = node.label.size,
				edge.label = edge.labels,
				edge.label.cex = edge.label.size,
				vertex.color = V(igraph)$color, 
				edge.color = E(igraph)$color,
				edge.width = edge.size, 
				vertex.size = node.size,
				edge.arrow.mode = arrow.shape,
				layout = layout_with_fr)
	} else {
		plot(igraph, edge.curve = 0.1, 
				vertex.label = node.labels,
				vertex.label.cex = node.label.size,
				vertex.color = V(igraph)$color, 
				edge.color = E(igraph)$color, 
				edge.width = edge.size, 
				vertex.size = node.size, 
				layout = layout_with_fr)
	}
	
		
	legend(position.legend.node, legend = names(color.nodes), fill = color.nodes, title = col.nodes)
	legend(position.legend.edge, legend = names(color.edges), fill = color.edges, title = col.edges)
	return
}