library(ggplot2)
library(glue)
library(dplyr)
library(RColorBrewer)

source('functions.network.R')

idh = "mut"
cohort <- "TCGA"
dims = "10"

input.file <- glue("ICA-JADE-{dims}.interactions.lr.filter-IDH{idh}-tcga.cgga-technical.txt")
output.file <- glue("ICA-JADE-{dims}.interactions.lr.filter.simultaneous.Status-IDH{idh}-tcga.cgga-technical.txt")
interactions <- read.table(input.file, sep = "\t", header = T)
interactions$type <- ifelse(interactions$prognosis == "Better", "Anti-tumor", "Pro-tumor")



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

bins <- interactions$bin
direction.genes <- interactions$direction.genes

simultaneous.state.lrs <- get.simultaneous.status(bins, direction.genes)
interactions$state.lr <- simultaneous.state.lrs
interactions$lr.activity <- ifelse(interactions$state.lr == "Both Expressed", "Active", "Inactive")

interactions.lrs.unique <- unique(interactions[ ,c('partners', 'type', 'lr.activity')])
table(interactions.lrs.unique$type, interactions.lrs.unique$lr.activity)

table(interactions$type, interactions$lr.activity)

write.table(interactions, file = output.file, sep = "\t", quote = F, row.names = F)


