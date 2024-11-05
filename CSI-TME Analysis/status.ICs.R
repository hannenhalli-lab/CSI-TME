library(ggplot2)
library(glue)
library(dplyr)
library(RColorBrewer)


idh = "mut"
cohort <- "TCGA"
dims = "10"




get.simultaneous.status <- function(bins, direction.genes) {
	####assign increase or decrease status to the interacting partners based on direction in IC and their interaction bin
	
	code.direction <- c("positive" = 1, "negative" = -1)
	code.bins <- c("Bin.1" = -1, "Bin.3" = 0, "Bin.9" = 1)
	
	f <- function(direction.genes) {
		###this will assign the joint initial status based on the two direction of IC###
		direction.genes.list <- strsplit(direction.genes, ":")
		domain.f <- lapply(direction.genes.list, function(x) {c1 = code.direction[x[[1]]]; c2 = code.direction[x[[2]]]; c(c1, c2)})
		change.initial <- lapply(domain.f, function(x) sum(x)) %>% unlist()
		return(change.initial)
	}
	
	g <- function(change.initial, bins) {
		###this will modify the joint initial status based on the bin of the interaction###
		bin.values <- code.bins[bins]
		change.final <- change.initial * bin.values
		return(change.final)
	}
	
	change.initial <- f(direction.genes)
	change.final <- g(change.initial, bins)
	
	###name the final status####
	change.final <- ifelse(change.final == 2, "Both Expressed", ifelse(change.final == -2, "Both Knocked", "One Expressed"))
	###additional correction for the cased in Bin.3#####
	change.final <- ifelse(bins == "Bin.3", ifelse(direction.genes == "positive:negative", "Both Expressed", "One Expressed"), change.final)
	
	return(change.final)
	
}