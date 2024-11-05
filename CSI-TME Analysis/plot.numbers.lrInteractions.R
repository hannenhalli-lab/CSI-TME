library(ggplot2)
library(glue)
library(dplyr)
library(pheatmap)

lr.file <- "ICA-JADE-10.interactions.lr.filter.simultaneous.Status-IDHmut-tcga.cgga-technical.txt"
output.file <- "plots/heatmap.lr.active.inactive.svg"
interactions <- read.table(lr.file, sep = "\t", header = T)
interactions$type <- ifelse(interactions$prognosis == "Better", "Anti-tumor", "Pro-tumor")

#data.network <- interactions[ ,c("lr.part1", "lr.part2", "bin", "type")]
#names(data.network) <- c("partner1", "partner2", "bin", "type")

data.network <- interactions[ ,c("lr.part1", "lr.part2", "lr.activity", "type")]
names(data.network) <- c("partner1", "partner2", "lr.activity", "type")


data.plot <- table(data.network$lr.activity, data.network$type)

svglite::svglite(file = output.file, height = 3, width = 3)
pheatmap(data.plot, display_numbers = T, fontsize_number = 12, cluster_row = F, cluster_col = F)
dev.off()  


