library(ggplot2); library(ggpubr)

input.file <- "mutationLoad.InteractionLoad.TCGA.txt"

load.data <- read.table(input.file, sep = "\t", header = T)

data.plot <- reshape2::melt(load.data, id.vars = c("Row.names", "x"))





data.plot$variable <- gsub("\\.", "-", data.plot$variable)
data.plot$variable <- factor(data.plot$variable, levels = c("Pro-Tumor", "Anti-Tumor", "total"))
scatter.plot <- ggplot(data = data.plot, aes(x = value, y = x)) + geom_point() + geom_smooth(method = "lm") + 
						ggpubr::stat_cor(method = "spearman", col = "red") + facet_wrap(~variable, scale = "free") + theme_bw() +
						xlab("# Interactions") + ylab("# Mutated genes")
						
						
output.file = "plots/ICA-JADE-10.FDR0.20.robustInteractions-IDHmut-mutationLoad.interactionLoad.svg"
svglite::svglite(file = output.file, height = 3, width = 9)
scatter.plot
dev.off()



data.plot <- data.plot[data.plot$variable == "total", ]


scatter.plot <- ggplot(data = data.plot, aes(x = value, y = x)) + geom_point() + geom_smooth(method = "lm") + 
						ggpubr::stat_cor(method = "spearman", col = "red") + theme_bw() +
						xlab("# Interactions") + ylab("# Mutated genes")
						
						
output.file = "plots/ICA-JADE-10.FDR0.20.robustInteractions-IDHmut-mutationLoad.interactionLoad-total.svg"
svglite::svglite(file = output.file, height = 3, width = 3)
scatter.plot
dev.off()

