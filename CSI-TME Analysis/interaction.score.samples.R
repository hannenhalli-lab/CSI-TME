

########load the ICA projections from the GLASS cohort##########
########call the activity ICs bins of the GLASS cohort##########
########call the interactions bins of the GLASS cohort##########
########for each interaction, observe the change wrt pre- or post- condition##########


#library(clusterProfiler)
library(dplyr); library(magrittr);
library(reshape2); library(ggplot2);
library(glue); 

source('../functionsForInteractionAnalysis.R')




idh <- "mut"; dims = 10; cohort = "GLASS"

input.name <- glue("ICA-JADE-{dims}.cellTypes.{tolower(cohort)}.IDH{idh}-filtered-noNorm-optimalRank.Rda")
load(input.name)

binMapList <- componentBinningFunction(pcaList = ica.res.list, factorization = "ICA")


interaction.file <- glue("ICA-JADE-{dims}.FDR0.25.robustInteractions-IDH{idh}-tcga.cgga-technical.txt")
interaction.data <- read.table(interaction.file, sep = "\t", header = T)
names(interaction.data)
interaction.data$prognosis <- ifelse(interaction.data$HR.x > 0, "Worse", "Better")
interactions.robust <- interaction.data %>% subset(., FDR.x < 0.20 & Reproduciblity > 7 & prognosis.x == prognosis.y & pvalue.y < 0.05)

progressIndex = 0
sample.matrix.list <- list()
for (index in 1:nrow(interactions.robust)) {
	interaction.string <- interactions.robust[index, ]
	
	cellTypes <- strsplit(interaction.string$cellTypes.x, ":") %>% unlist()
	components <- strsplit(interaction.string$components.x, ":") %>% unlist()
	int.bin <- interaction.string$bin.x
	
	cellType1 <- cellTypes[1]
	cellType2 <- cellTypes[2]
	
	dim1 <- components[1]
	dim2 <- components[2]
	
	comp1.bins <- binMapList[[cellType1]][ ,dim1]
	comp2.bins <- binMapList[[cellType2]][ ,dim2]
	
	binValues <- data.frame(barcodes = names(comp1.bins), bin1 = comp1.bins, bin2 = comp2.bins) %>% set_rownames(NULL)
	
	if (int.bin == "Bin.1") {
		bins <- ifelse(binValues$bin1 == 0 & binValues$bin2 == 0, 1, 0)
	}
	if (int.bin == "Bin.3") {
		bins <- ifelse(binValues$bin1 == 2 & binValues$bin2 == 0, 1, 0)
	}
	if (int.bin == "Bin.9") {
		bins <- ifelse(binValues$bin1 == 2 & binValues$bin2 == 2, 1, 0)
	}
	sample.matrix.list[[index]] <- bins %>% set_names(names(comp1.bins)) %>% data.frame()
}


sample.matrix.df <- do.call("cbind", sample.matrix.list) %>% set_colnames(interactions.robust$string) %>% t()

sample.matrix.df <- sample.matrix.df %>% data.matrix %>% reshape2::melt()
sample.matrix.df$Var2 <- gsub("\\.", "-", sample.matrix.df$Var2)
colnames(sample.matrix.df) <- c("string", "samples", "value")
sample.matrix.df$samples <- gsub("\\.", "-", sample.matrix.df$samples) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)

#write.table(sample.matrix.df, glue("interaction.samples.distribution.IDH{idh}-{cohort}.txt"), sep = "\t", quote = F)

clinical.file <- "glass_metadata_combined.csv"
clinical <- read.table(clinical.file, sep = ",", header = T)
clinical <- clinical %>% subset(., select = c(sample = sample_barcode, Tumor_Type, idh_status))
colnames(clinical) <- c("samples", "type", "idh")





######plots##########
###################fisher like##########
data.plot <- merge(sample.matrix.df, clinical, by = "samples")
data.plot$samples <- gsub("-\\w+$", "", data.plot$samples)
data.plot <- merge(data.plot, interactions.robust, by = "string")
data.plot$int.type <- ifelse(data.plot$prognosis == "Worse", "Pro-tumor", "Anti-tumor")


data.list <- split(data.plot, data.plot$string)

fisher.res.list <- list(); index = 0
for (interaction in names(data.list)) {
	index <- index + 1
	data.interaction <- data.list[[interaction]]
	int.type <- data.interaction$int.type[1]
	mat <- table(data.interaction$value, data.interaction$type)
	if (nrow(mat) == 2 & ncol(mat) == 2) {
		mat <- mat[c(2, 1), c(2, 1)]
		fisher.res <- fisher.test(mat)
		estimate <- fisher.res$estimate; pvalue <- fisher.res$p.value
	} else {
		estimate <- NA; pvalue <- NA
	}
	
	fisher.res.list[[interaction]] <- data.frame(interaction = interaction, oddsRatio = estimate, pvalue = pvalue, prognosis = int.type)
}
fisher.res.df <- bind_rows(fisher.res.list, .id = "string")
fisher.res.df$FDR <- p.adjust(fisher.res.df$pvalue, method = "fdr")

fisher.res.df <- fisher.res.df[fisher.res.df$pvalue < 0.05, ] %>% na.omit()
fisher.res.df$oddsRatio <- ifelse(is.finite(fisher.res.df$oddsRatio), fisher.res.df$oddsRatio, 7)
high.sig <- fisher.res.df$string



################penetrance###############

#data.plot <- merge(sample.matrix.df, clinical, by = "samples")
#data.plot$samples <- gsub("-\\w+$", "", data.plot$samples)

#data.plot <- merge(data.plot, interactions.robust, by = "string")

data.plot <- data.plot[ ,c('string', 'value', 'type', 'int.type')]

data.plot <- data.plot %>% group_by(string, type, int.type) %>% summarize(proportion = sum(value)/length(value))

box.plot <-  ggplot(data = data.plot, aes(x = type, y = proportion)) + geom_boxplot() + facet_wrap(~int.type) + theme_bw() + ggpubr::stat_compare_means()

data.plot <- data.plot %>% reshape2::dcast(string + int.type ~ type, value.var = 'proportion', fun.aggregate = mean, na.rm = T, na.action = "na.pass") ## fun.aggregate operates to average multiple recurrent samples



data.plot$remodelling <- ifelse(data.plot$string %in% high.sig, "P < 0.05", "P > 0.05")


scatter.plot <- ggplot(data = data.plot, aes(x = Primary, y = Recurrent, col = remodelling)) + geom_point() + facet_wrap(~int.type) + theme_bw() + 
		coord_cartesian(xlim = c(0,0.30), ylim = c(0,0.30)) + geom_abline(slope = 1, coef = 0)
#data.heatmap <- data.plot %>% reshape2::dcast(string + int.type ~ type, value.var = 'proportion', fun.aggregate = mean, na.rm = T, na.action = "na.pass") ## fun.aggregate operates to a




output.file <- glue("plots/scatterplot.primary.recurrent.penetrance.{cohort}.IDH{idh}.svg")
svglite::svglite(file = output.file, width = 5, height = 2.5)
scatter.plot
dev.off()


output.file <- glue("plots/boxplot.primary.recurrent.penetrance.{cohort}.IDH{idh}.svg")
svglite::svglite(file = output.file, width = 5, height = 2.5)
box.plot
dev.off()



####################interactionload##################
data.plot <- merge(sample.matrix.df, clinical, by = "samples")
data.plot$samples <- gsub("-\\w+$", "", data.plot$samples)
data.plot <-  merge(data.plot, interactions.robust, by = "string")
data.plot$int.type <- ifelse(data.plot$prognosis == "Worse", "Pro-tumor", "Anti-tumor")


primary.samples <- data.plot[data.plot$type == "Primary", 'samples'] %>% unique()
recurrent.samples <- data.plot[data.plot$type == "Recurrent", 'samples'] %>% unique()
common <- intersect(primary.samples, recurrent.samples)



#data.plot <- data.plot[data.plot$samples %in% common, ]

data.plot <- data.plot %>% group_by(samples, type, int.type) %>% summarize(proportion = sum(value)/length(value))


box.plot <- ggplot(data = data.plot, aes(x = type, y = proportion)) + geom_boxplot() + facet_wrap(~int.type) + theme_bw() + ggpubr::stat_compare_means()

data.plot <- data.plot %>% reshape2::dcast(samples + int.type ~ type, value.var = 'proportion', fun.aggregate = mean, na.rm = T, na.action = "na.pass") ## fun.aggregate operates to average multiple recurrent samples


scatter.plot <- ggplot(data = data.plot, aes(x = Primary, y = Recurrent)) + geom_point() + facet_wrap(~int.type) + theme_bw() +
	coord_cartesian(xlim = c(0, 0.55), ylim = c(0,0.55)) + geom_abline(slope = 1, intercept = 0) 
	
	
	
#data.plot.better <- data.plot[data.plot$prognosis == "Better", ]
#data.plot.worse <- data.plot[data.plot$prognosis == "Worse", ]


output.file <- glue("plots/scatterplot.primary.recurrent.interactionLoad.{cohort}.IDH{idh}.svg")
svglite::svglite(file = output.file, width = 5, height = 2.5)
scatter.plot
dev.off()


output.file <- glue("plots/boxplot.primary.recurrent.interactionLoad.{cohort}.IDH{idh}.svg")
svglite::svglite(file = output.file, width = 5, height = 2.5)
box.plot
dev.off()



















##########################obsolete##########################


#########plot results#########

cohort <- "GLASS"
sample.matrix.df <- read.table(glue("interaction.samples.distribution.IDH{idh}-{cohort}.txt"), sep = "\t", row.names = 1)

sample.matrix.df <- sample.matrix.df %>% data.matrix %>% reshape2::melt()
sample.matrix.df$Var2 <- gsub("\\.", "-", sample.matrix.df$Var2)
colnames(sample.matrix.df) <- c("string", "samples", "value")
sample.matrix.df$samples <- gsub("\\.", "-", sample.matrix.df$samples) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)



if (cohort == "TCGA") {
	cancerType <- c("GBM", "LGG")
	clinical <- read.table("AllAvialableClinicalDataMod.txt", sep = "\t", header = T)
	clinical <- clinical[clinical$type %in% cancerType, ]
	
	codel <- read.csv("../All_metadata-CCGA.csv", header = T) %>% set_names(c('samples' ,'project', 'time', 'status', 'IDH', 'grade', 'age', 'subtype', 'codel'))
	codel <- codel[grep("TCGA", codel$samples), c('IDH', 'samples', 'codel', 'grade')]
	codel$samples <- gsub("\\.", "-", codel$samples) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
	clinical <- merge(clinical, codel, by = "samples")
}

if (cohort == "CGGA") {
	clinical <- read.csv("../All_metadata-CCGA.csv", header = T) %>% set_names(c('samples' ,'project', 'time', 'status', 'IDH', 'grade', 'age', 'subtype', 'codel'))
	gender <- read.table("cgga.combinedMetadata.txt", sep = "\t", header = T) %>% .[ ,c('samples', 'Gender')]
	clinical <- merge(clinical, gender, by = "samples")
	names(clinical)[names(clinical) == 'Gender'] <- "sex"
}

if (cohort == "GLASS") {
	clinical.file <- "glass_metadata_combined.csv"
	clinical <- read.table(clinical.file, sep = ",", header = T)
	clinical <- clinical %>% subset(., select = c(sample = sample_barcode, Tumor_Type, idh_status))
	colnames(clinical) <- c("samples", "type", "idh")
}


data.plot <- merge(sample.matrix.df, clinical, by = "samples")
data.plot$samples <- gsub("-\\w+$", "", data.plot$samples)

primary.samples <- data.plot[data.plot$type == "Primary", 'samples'] %>% unique()
recurrent.samples <- data.plot[data.plot$type == "Recurrent", 'samples'] %>% unique()
common <- intersect(primary.samples, recurrent.samples)



###################fisher like##########
data.plot <- merge(sample.matrix.df, clinical, by = "samples")
data.plot$samples <- gsub("-\\w+$", "", data.plot$samples)
#data.plot <- data.plot[data.plot$samples %in% common, ]

data.plot <- merge(data.plot, interactions.robust, by = "string")
data.list <- split(data.plot, data.plot$string)



fisher.res.list <- list(); index = 0
for (interaction in names(data.list)) {
	index <- index + 1
	data.interaction <- data.list[[interaction]]
	prognosis <- data.interaction$prognosis[1]
	mat <- table(data.interaction$value, data.interaction$type)
	if (nrow(mat) == 2 & ncol(mat) == 2) {
		mat <- mat[c(2, 1), c(2, 1)]
		fisher.res <- fisher.test(mat)
		estimate <- fisher.res$estimate; pvalue <- fisher.res$p.value
	} else {
		estimate <- NA; pvalue <- NA
	}
	
	fisher.res.list[[interaction]] <- data.frame(interaction = interaction, oddsRatio = estimate, pvalue = pvalue, prognosis = prognosis)
}
fisher.res.df <- bind_rows(fisher.res.list, .id = "string")
fisher.res.df$FDR <- p.adjust(fisher.res.df$pvalue, method = "fdr")

fisher.res.df <- fisher.res.df[fisher.res.df$pvalue < 0.05, ] %>% na.omit()
fisher.res.df$oddsRatio <- ifelse(is.finite(fisher.res.df$oddsRatio), fisher.res.df$oddsRatio, 7)
high.sig <- fisher.res.df$string




################paired interaction###############
data.plot <- merge(sample.matrix.df, clinical, by = "samples")
data.plot$samples <- gsub("-\\w+$", "", data.plot$samples)

data.plot <- merge(data.plot, interactions.robust, by = "string")
data.plot$int.type <- ifelse(data.plot$prognosis == "Worse", "Pro-tumor", "Anti-tumor")

data.plot <- data.plot[ ,c('string', 'value', 'type', 'int.type')]

data.plot <- data.plot %>% group_by(string, type, int.type) %>% summarize(proportion = sum(value)/length(value))

data.heatmap <- data.plot %>% reshape2::dcast(string + int.type ~ type, value.var = 'proportion', fun.aggregate = mean, na.rm = T, na.action = "na.pass") ## fun.aggregate operates to average multiple recurrent samples

data.heatmap.plot <- data.heatmap[ ,c(3:4)] %>% set_rownames(data.heatmap$string) %>% data.matrix()
data.rowannotation <- data.heatmap %>% subset(., select = int.type) %>% set_rownames(data.heatmap$string) #%>% data.matrix()
data.rowannotation$pvalue <- ifelse(rownames(data.heatmap.plot) %in% high.sig, "<0.05", "other")


pheat.output <- glue("plots/pheatamp.interactions.proportion.{cohort}.IDH{idh}.svg")
library(gridExtra)
#####following three lines from https://gist.github.com/jakalssj3/3f24e0a8e18f0badb4f53de487474dc5######
my_pheatmap <- pheatmap::pheatmap(data.heatmap.plot, cluster_cols = F, annotation_row = data.rowannotation, clustering_method = "ward.D2", fontsize_row = 5)[[4]]
my_plot <- grid.arrange(my_pheatmap, nrow = 1, ncol = 1)
ggsave(filename = pheat.output, plot = my_plot, height = 6, width = 3)
dev.off()

######################################################################################
######################################################################################
####################################paired samples####################################
data.plot <- merge(sample.matrix.df, clinical, by = "samples")
data.plot$samples <- gsub("-\\w+$", "", data.plot$samples)

data.plot <-  merge(data.plot, interactions.robust, by = "string")

primary.samples <- data.plot[data.plot$type == "Primary", 'samples'] %>% unique()
recurrent.samples <- data.plot[data.plot$type == "Recurrent", 'samples'] %>% unique()
common <- intersect(primary.samples, recurrent.samples)



data.plot <- data.plot[data.plot$samples %in% common, ]

data.plot <- data.plot %>% group_by(samples, type, prognosis) %>% summarize(proportion = sum(value)/length(value))

data.plot <- data.plot %>% reshape2::dcast(samples + prognosis ~ type, value.var = 'proportion', fun.aggregate = mean, na.rm = T, na.action = "na.pass") ## fun.aggregate operates to average multiple recurrent samples


data.plot.better <- data.plot[data.plot$prognosis == "Better", ]
data.plot.worse <- data.plot[data.plot$prognosis == "Worse", ]


output.file <- glue("plots/scatterplot.primary.recurrent.interactionLoad-Protumor.{cohort}.IDH{idh}.svg")
svg(file = output.file, width = 3, height = 3)
ggplot(data = data.plot.worse, aes(x = Primary, y = Recurrent)) + geom_point() + coord_cartesian(xlim = c(0, 0.5), ylim = c(0,0.5)) + geom_abline(slope = 1, intercept = 0) + theme_bw()
dev.off()

output.file <- glue("plots/scatterplot.primary.recurrent.interactionLoad-Antitumor.{cohort}.IDH{idh}.svg")
svg(file = output.file, width = 3, height = 3)
ggplot(data = data.plot.better, aes(x = Primary, y = Recurrent)) + geom_point() + coord_cartesian(xlim = c(0, 0.5), ylim = c(0,0.5)) + geom_abline(slope = 1, intercept = 0) + theme_bw()
dev.off()




data.plot.worse$delta <- data.plot.worse$Recurrent - data.plot.worse$Primary
data.plot.worse$direction <- ifelse(data.plot.worse$delta >= 0.05, "Increased", ifelse(data.plot.worse$delta <= -0.05, "Decreased", "NoChange"))


data.plot.better$delta <- data.plot.better$Recurrent - data.plot.better$Primary
data.plot.better$direction <- ifelse(data.plot.better$delta >= 0.05,  "Increased", ifelse(data.plot.better$delta <= -0.05, "Decreased", "NoChange"))
