library(dplyr); library(glue); library(magrittr);
library(ggplot2); library(ggpubr); library(data.table)

cellTypes <- c("Bcell", "Endothelial", "Malignant", "Myeloid", "Oligos", "Stromal", "Tcell")

direction <- "positive"
for (cellType in cellTypes) {
	input <- glue("singleCellDistances/distances_PC_singleCell_idhMut-{cellType}-{direction}Genes.Rda")
	load(input)
	
	data.plot <- rbindlist(distances_list)
	data.plot <- data.plot[data.plot$Label != "Within Non-members", ]
	p <- ggplot(data = data.plot, aes(x = Label, y = Distance)) + geom_boxplot() + theme_bw() + 
			facet_wrap(~IC, nrow = 2) + stat_compare_means(label = "p.format") + 
			theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1))
			
	outputFile <- glue("distances_PC_singleCell_cellType}-{direction}Genes-plot.svg")
	svglite::svglite(file = outputFile, width = 6, height = 4)
	p
	dev.off()
}

#####median distances#####

distance_pvalues_list <- list()
direction <- "positive"
for (cellType in cellTypes) {
	input <- glue("singleCellDistances/distances_PC_singleCell_idhMut-{cellType}-{direction}Genes.Rda")
	message(glue("loading data for {cellType}\n"))
	load(input)
	
	data.plot <- rbindlist(distances_list)
	data.plot <- data.plot[data.plot$Label != "Within Non-members", ]
	distances <- data.plot %>% group_by(Label, IC) %>% summarize(distance = median(Distance, na.rm = T))
	
	message(glue("computing p-values for distances\n"))
	stat.test <- compare_means(
	 	Distance ~ Label, 
		group.by = "IC",
		data = data.plot,
	 	method = "wilcox.test"
	)
	
	distances <- reshape2::dcast(distances, IC ~ Label, value.var = "distance")
	colnames(distances) <- gsub(" .+", "", colnames(distances))
	
	distances <- merge(distances, stat.test, by = "IC")
	distances$percent <- ((distances$Between - distances$Within) / distances$Between) * 100
	
	distance_pvalues_list[[cellType]] <- distances
	message(glue("done!\n"))
}

outputFile <- glue("distances_PC_singleCell-{direction}Genes-pvalues.txt")
distance_pvals <- bind_rows(distance_pvalues_list, .id = "cellType")
write.table(distance_pvals, outputFile, sep = "\t", quote = F, row.names = F)



######makeplot#####
library(glue); library(dplyr); library(ggplot2); 
library(circlize); library(magrittr)

add_categories <- function(df, group_column, category_column, categories) {
  # Create all combinations of the group and category
  all_combinations <- expand.grid(
     group = unique(df[[group_column]]),
     category = categories,
     stringsAsFactors = FALSE
   )
   names(all_combinations) <- c(group_column, category_column)
  
   by_join <- setNames(nm = c(group_column, category_column)) %>% set_names(NULL)
   names(by_join) <- NULL
  
   # Left join the original data frame with all combinations
   # This ensures all possible combinations are present, filling missing ones with NA
   #result <- all_combinations %>% left_join(df, by = ≈)
   
   result <- merge(df, all_combinations, by = by_join, all = T)   
  
  return(result)
}


direction <- "positive"
inputFile <- glue("distances_PC_singleCell-{direction}Genes-pvalues.txt")


input.data.positive <- read.table(inputFile, sep = "\t", header = T)
input.data.positive$IC <- gsub("Dim", "IC", input.data.positive$IC)
all.ICs <- paste("IC.", 1:10, sep = "")

input.data.positive <- add_categories(df = input.data.positive, group_column = "cellType", category_column = "IC", categories = all.ICs)
input.data.positive$colors <- ifelse(input.data.positive$percent > 10 & input.data.positive$p.adj < 0.20, "darkgreen", "darkred")



direction <- "negative"
inputFile <- glue("distances_PC_singleCell-{direction}Genes-pvalues.txt")


input.data.negative <- read.table(inputFile, sep = "\t", header = T)
input.data.negative$IC <- gsub("Dim", "IC", input.data.negative$IC)
all.ICs <- paste("IC.", 1:10, sep = "")
cellTypes <- unique(input.data.negative$cellType)

input.data.negative <- add_categories(df = input.data.negative, group_column = "cellType", category_column = "IC", categories = all.ICs)
input.data.negative$colors <- ifelse(input.data.negative$percent > 10 & input.data.negative$p.adj < 0.20, "darkgreen", "darkred")


xlim.matrix <- table(input.data.positive$cellType) %>% as.matrix()
xlim.matrix <- cbind(0, xlim.matrix)


#data.within <- split(input.data.positive, input.data.positive$cellType) %>% lapply(., function(x) x$Within)
#data.between <- split(input.data.positive, input.data.positive$cellType) %>% lapply(., function(x) -x$Between) ###negative to plot on inner side


data.percent.positive <- split(input.data.positive, input.data.positive$cellType) %>% lapply(., function(x) x$percent)
data.percent.negative <- split(input.data.negative, input.data.negative$cellType) %>% lapply(., function(x) x$percent)

color.positive <- split(input.data.positive, input.data.positive$cellType) %>% lapply(., function(x) x$colors)
color.negative <- split(input.data.negative, input.data.negative$cellType) %>% lapply(., function(x) x$colors)

position.label.x <- seq(1, by = 1, length.out = length(cellTypes))


output.file <- "pc_distances_plot.svg"

svglite::svglite(output.file, width = 7, height = 7)

circos.par(cell.padding=c(0,0,0,0), start.degree = 90, gap.degree = 6)
circos.initialize(sectors = cellTypes, xlim = xlim.matrix)
#ylim <- range(unlist(data.percent.positive), na.rm = T)

ylim <- lapply(data.percent.positive, range, na.rm = T) %>% do.call("rbind", .)
ylim[ ,2] <- ylim[ ,2] + 5

max.y <- max(ylim) + 15
mean.x <- length(all.ICs) / 2
position.label.x <- rep(mean.x, length(cellTypes))
#position.label.y <- rep(max.y, length(cellTypes))
position.label.y <- ylim[ ,2] + 4
position.label.matrix <- cbind(position.label.x, position.label.y) %>% set_rownames(cellTypes)


circos.track(ylim = ylim, y = 1:length(cellTypes), panel.fun = function(x, y) {
	value <- data.percent.positive[[y]]
	color.vector <- color.positive[[y]]
	name = get.cell.meta.data("sector.index")
	x.mat <- c(1, 1) %>% matrix(ncol = 2)
	xaxis.position <- seq.int(from = 0.5, to = 9.5, length.out = 10)
	xaxis.label <- c(1,10,2:9)
	sector.label.position <- position.label.matrix[y, ] %>% matrix(., ncol = 2)
	circos.barplot(value, pos = xaxis.position, col = color.vector, lwd = 0.05, lty = 2)
	circos.yaxis(side = "right", labels.cex = 0.7)	
	circos.axis(sector.index = name, h = "bottom", labels = xaxis.label, labels.niceFacing = T, direction = "inside", labels.cex = 0.85, major.at = xaxis.position)
	circos.text(labels = name,  x = sector.label.position, facing = "outside", niceFacing = TRUE, cex = 2, adj = c(0.5, 0.5))
	
})

set_track_gap(gap = 0.15)

#ylim <- range(unlist(data.percent.negative), na.rm = T)
ylim <- lapply(data.percent.negative, range, na.rm = T) %>% do.call("rbind", .)
ylim[ ,2] <- ylim[ ,2] + 5
circos.track(ylim = ylim, y = 1:length(cellTypes), panel.fun = function(x, y) {
	value <- data.percent.negative[[y]]
	color.vector <- color.negative[[y]]
	xaxis.position <- seq.int(from = 0.5, to = 9.5, length.out = 10)
	xaxis.label <- c(1,10,2:9)
	circos.barplot(value, pos = xaxis.position, col = color.vector)
	
	circos.yaxis(side = "right", labels.cex = 0.7)
	#circos.xaxis()
})
dev.off()
####following also works
#circos.track(ylim = ylim, factors = input.data.positive$cellType, x = input.data.positive$IC, y = input.data.positive$percent,
#             panel.fun = function(x, y) {
#               circos.barplot(y, 1:10 - 0.5, col = 1:10)
#             })
			 
