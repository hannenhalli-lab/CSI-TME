###########plot#############
library(dplyr); library(magrittr); library(glue)
library(pheatmap); library(ggplot2); library(tidyverse)
source('../functionsForInteractionAnalysis.R')


input.type <- "interactions"
dims <- 10; idh = "mut"
if (input.type == "interactions") {
	file.input <- glue("ICA-JADE-{dims}.FDR0.20.robustInteractions-IDH{idh}-association.Mutations.txt")
}
if (input.type == "cellTypes") {
	file.input <- glue("ICA-{method}-{dims}.cellTypes.{cohort}.IDH{idh}-assocaition.Mutations.txt")
}

create.colors <- function(min.bound = 0, max.bound = 10, middle.value = 1, partition.length = 10) {
	#min.bound --> minimum value in the data to color for
	#max.bound --> maximum value in the data to color for
	#middle.value --> middle value to center around the colors
	#partition.length <- how many intervals each values should be partininted into below and aboce the middle.value
	
	
	by.lower <- middle.value / partition.length
	by.upper <- max.bound / partition.length

	breaks.lower <- seq(min.bound, middle.value, by = by.lower)
	breaks.upper <- seq(middle.value, max.bound, by = by.upper)[-1]

	breaks <- c(breaks.lower, breaks.upper)
	paletteLength <- length(breaks) - 1
	color.pallete <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
	return(list(colors = color.pallete, breaks = breaks))
	
}


input.data <- read.table(file.input, sep = "\t", header = T)
input.data$Association <- ifelse(input.data$Estimate > 1 , "Positive", "Negative") ##### estimate is odds ratio

association.significant <- input.data[input.data$FDR < 0.20,  ]

data.plot <- reshape2::dcast(data = association.significant, formula = predictor ~ response, value.var = "Estimate") %>% column_to_rownames(var = "predictor")
data.plot[is.na(data.plot)] <- 1 ####no association means odds ratio == 1

allValues <- unlist(data.plot)

data.plot[data.plot == -Inf] <- min(allValues[is.finite(allValues)], na.rm = TRUE)
data.plot[data.plot == +Inf] <- max(allValues[is.finite(allValues)], na.rm = TRUE)

rowAnnotationDf <- association.significant %>% subset(., select = c(predictor, prognosis)) %>% unique() %>% set_rownames(.$predictor)
rowAnnotationDf$int.type <- ifelse(rowAnnotationDf$prognosis == "Better", "Anti-tumor", "Pro-tumor")
rowAnnotationDf <- rowAnnotationDf %>% subset(., select = c(int.type))


outputFile <- glue("plots/pheatmap.mutations.interactions.associations.svg")

my.pallete <- create.colors(min.bound = 0, max.bound = 20, middle.value = 1, partition.length = 10)

plot.pheatmap <- pheatmap(data.plot, annotation_row = rowAnnotationDf, color = my.pallete[['colors']], breaks = my.pallete[['breaks']], fontsize_row = 8, fontsize_col = 8)

svglite::svglite(outputFile, width = 7, height = 12)
plot.pheatmap
dev.off()
