library(TCGAbiolinks); library(glue); 
library(dplyr); library(data.table);
library(maftools);  library(magrittr)
source('../functionsForInteractionAnalysis.R')

###
##options(pkgType="binary")
##install.packages("SMDIC")

association.function <- function(data1, data2, covariates = NULL, method.regression) {
	library(simsalapar)
	message("currently, how to deal with NA values in not added in this function")
	list.variable1 <- colnames(data1)
	list.variable2 <- colnames(data2)
	variable.combinations <- expand.grid(list.variable1, list.variable2)
	size.iteration <- nrow(variable.combinations)
	
	name.output <- paste(variable.combinations[ ,1], variable.combinations[ ,2], sep = "-")
	coefficient.list <- list()
	for (index in 1:size.iteration) {
		variable1 <- variable.combinations[index, 1]
		variable2 <- variable.combinations[index, 2]
		
		test.data1 <- data1[ ,variable1] %>% data.frame() %>% set_colnames("response") %>% mutate(samples = rownames(data1))
		test.data2 <- data2[ ,variable2] %>% data.frame() %>% set_colnames("predictor") %>% mutate(samples = rownames(data2))
		
		if (!is.null(covariates)) {
			input.list <- list(test.data1, test.data2, covariates)
		} else {
			input.list <- list(test.data1, test.data2)
		}
		data.test <- input.list %>% purrr::reduce(full_join, by='samples')
		data.test <- subset(data.test, select = -samples)
		
		if (method.regression == "fisher") {
			data.test <- na.omit(data.test) ###remove where NA are present###
			regression.output <-  tryCatch.W.E(table(data.test) %>% fisher.test)
			coefficents <- data.frame(Estimate = regression.output$value$estimate, p.value = regression.output$value$p.value)
		} else {
			regression.output <- tryCatch.W.E(glm(response ~ predictor + ., data = data.test, family = method.regression))
			#regression.output <- glm(response ~ predictor + ., data = data.test, family = method.regression)
			coefficents <- summary(regression.output$value)[['coefficients']] %>% data.frame
			colnames(coefficents) <- c("Estimate", "StdError", "t.value", "p.value")
		}
		warning <- regression.output$warning
		
		if (is.null(warning)) {
			warning <- "No"
		} else {
			warning <- "Yes"
		}
		
		coefficents$warning <- warning
		coefficents$response <- variable1
		coefficents$predictor <- variable2
		
		if (method.regression == "fisher") {
			coefficient.list[[index]] <- coefficents
		} else {
			coefficient.list[[index]] <- coefficents[-1, ] ####first row is intercept
		}
	}
	names(coefficient.list) <- name.output
	return(coefficient.list)
}

maf.input <- "glioma.mutation_maftools.maf"

maf.file.data <- fread(maf.input, nThread = 1)
maf <- read.maf(maf.file.data)
mutation.data <- maf@data[ ,c('Hugo_Symbol', 'Tumor_Sample_Barcode')]
mutation.data <- dcast.data.table(mutation.data, Hugo_Symbol ~ Tumor_Sample_Barcode)
mutation.data <- data.frame(mutation.data)
rownames(mutation.data) <- mutation.data[ ,1]
mutation.data <- mutation.data[ ,-1]
colnames(mutation.data) <- gsub("\\.", "-", colnames(mutation.data)) %>% gsub("\\w-\\w+-\\w+-\\w+$", "", .)


idh = "mut"
cohort <- "tcga"
dims = 10; method = "JADE"
remove.bcells = F

####load patient clinical data#######
load(glue("ICA-{method}-{dims}.cellTypes.{cohort}.IDH{idh}-filtered-noNorm-optimalRank.Rda"))
interactionFile <- glue("ICA-JADE-{dims}.FDR0.25.robustInteractions-IDHmut-tcga.cgga-technical.txt")
interactions <- read.table(interactionFile, sep = "\t", header = T)
interactions.robust <- interactions %>% subset(., FDR.x < 0.20 & Reproduciblity > 7)


binMapList <- componentBinningFunction(pcaList = ica.res.list, factorization = "ICA", sample.names = NULL)
names(binMapList) <- gsub("\\.", "", names(binMapList))


######make unique sample names#########

for (cellType in names(binMapList)) {
	bin.map.data <- binMapList[[cellType]]
	sample.names <- gsub("\\w.\\w+.\\w+.\\w+$", "", rownames(bin.map.data) ) 
	sample.names <- make.names(sample.names, unique = TRUE) %>% gsub("\\.", "-", .)
	rownames(bin.map.data) <- sample.names
	binMapList[[cellType]] <- bin.map.data
}





if (remove.bcells == T) {
	indexes <- grep("Bcell", interactions.robust$string, invert = T)
	interactions.robust <- interactions.robust[indexes, ]
}


#interactions.robust <- interactions.robust[interactions.robust$prognosis.x == "Worst", ]
interaction.matrix.list <- list()
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
	
	interaction.matrix.list[[index]] <- bins %>% data.frame() %>% set_rownames(binValues$barcodes)

}
interaction.matrix <- do.call("cbind", interaction.matrix.list) 
colnames(interaction.matrix) <- interactions.robust$string

#####filter mutation data#####
sample.use <- rownames(interaction.matrix)
mutation.data.filtered <- mutation.data[ ,colnames(mutation.data) %in% sample.use] ###filter only IDH mutant
mutation.freqeuncy <- apply(mutation.data.filtered, 1, function(x) sum(x > 0) / length(x))
genes.use <- mutation.freqeuncy[mutation.freqeuncy > 0.01] %>% names
mutation.data.filtered <- mutation.data.filtered[genes.use, ] %>% t() %>% data.frame()
mutation.data.filtered[mutation.data.filtered > 0] <- 1


output.regression.list <- list()
for (cellType in names(binMapList)) {
	output.regression <- association.function(data1 = mutation.data.filtered, data2 = binMapList[[cellType]], method.regression = "binomial")
	output.regression <- bind_rows(output.regression, .id = "comparison")
	output.regression.list[[cellType]] <- output.regression
	message(glue("done for {cellType} cellType"))
}

output.regression <- bind_rows(output.regression.list, .id = "cellType")
output.regression$FDR <- p.adjust(output.regression$p.value, method = "fdr")
write.table(output.regression, file = glue("ICA-{method}-{dims}.cellTypes.{cohort}.IDH{idh}-assocaition.Mutations.txt"), sep = "\t", quote = F, row.names = F)



#output.regression <- association.function(data1 = mutation.data.filtered, data2 = interaction.matrix, method.regression = "binomial")
output.regression <- association.function(data1 = mutation.data.filtered, data2 = interaction.matrix, method.regression = "fisher")

output.regression <- bind_rows(output.regression, .id = "comparison")
output.regression$FDR <- p.adjust(output.regression$p.value, method = "fdr")

prognosis <- interactions.robust[ ,c('string', "prognosis.x")] %>% set_colnames(c('predictor', 'prognosis'))
output.regression <- merge(output.regression, prognosis, by = 'predictor')
output.regression <- merge(output.regression, mutation.freqeuncy, by.x = "response", by.y = 0)
colnames(output.regression)[9] <- "mutation.freqeuncy"
write.table(output.regression, file = glue("ICA-JADE-{dims}.FDR0.20.robustInteractions-IDH{idh}-association.Mutations.txt"), sep = "\t", quote = F, row.names = F)







###########################Mutation load###############################
library(ggplot2)
output.file <- glue("plots/ICA-JADE-{dims}.FDR0.20.robustInteractions-IDHmut-mutationLoad.interactionLoad.svg")
sample.use <- rownames(interaction.matrix)
mutation.data.filtered <- mutation.data[ ,colnames(mutation.data) %in% sample.use] ###filter only IDH mutant
mutation.data.filtered[mutation.data.filtered > 0] <- 1
mutation.load <- colSums(mutation.data.filtered)


#sample.use <- rownames(interaction.matrix)
#mutation.data.filtered <- mutation.data[ ,colnames(mutation.data) %in% sample.use] ###filter only IDH mutant
#mutation.freqeuncy <- apply(mutation.data.filtered, 1, function(x) sum(x > 0) / length(x))
#genes.use <- mutation.freqeuncy[mutation.freqeuncy > 0.01] %>% names
#mutation.data.filtered <- mutation.data.filtered[genes.use, ]
#mutation.data.filtered[mutation.data.filtered > 0] <- 1
#mutation.load <- colSums(mutation.data.filtered)


interactions.better <- interactions.robust[interactions.robust$prognosis.x == "Better", 'string']
interactions.worse <- interactions.robust[interactions.robust$prognosis.x == "Worst", 'string']


interaction.load.better <- rowSums(interaction.matrix[ ,colnames(interaction.matrix) %in% interactions.better])
interaction.load.worse <- rowSums(interaction.matrix[ ,colnames(interaction.matrix) %in% interactions.worse])

interaction.load <- data.frame(interaction.load.worse, interaction.load.better) %>% set_colnames(c("Pro.Tumor", "Anti.Tumor"))
interaction.load$total <- interaction.load$Pro.Tumor + interaction.load$Anti.Tumor

data.load <- merge(mutation.load, interaction.load, by = 0)
write.table(data.load, file = "mutationLoad.InteractionLoad.TCGA.txt", sep = "\t", quote = F, row.names = F)
data.load$Total <- data.load$Pro.Tumor + data.load$Anti.Tumor

data.plot <- reshape2::melt(data.load, id.vars = c("Row.names", "x"))
data.plot$variable <- gsub("load\\.", "", data.plot$variable)
data.plot$variable <- factor(data.plot$variable, levels = c("better", "worse", "total"))
scatter.plot <- ggplot(data = data.plot, aes(x = value, y = x)) + geom_point() + geom_smooth(method = "lm") + 
						ggpubr::stat_cor(method = "spearman", col = "red") + facet_wrap(~variable, scale = "free") + theme_bw() +
						xlab("# Interactions") + ylab("# Mutated genes")

svg(file = output.file, height = 3, width = 9)
scatter.plot
dev.off()



scatter.plot <- ggplot(data = data.plot, aes(x = value, y = x, col = variable)) + geom_point() + geom_smooth(method = "lm") + 
						ggpubr::stat_cor(method = "spearman", aes(color = variable)) +  theme_bw() +
						xlab("# Interactions") + ylab("# Mutated genes")






###########################interaction frequency###################################
output.file <- glue("plots/ICA-JADE-{dims}.FDR0.25.robustInteractions-IDHmut-frequency.prognosis.svg")
interaction.freq.better <- colSums(interaction.matrix[ ,colnames(interaction.matrix) %in% better]) / nrow(interaction.matrix)
interaction.freq.worse <- colSums(interaction.matrix[ ,colnames(interaction.matrix) %in% interactions.worse]) / nrow(interaction.matrix)

expectation <- nrow(interaction.matrix) / 9 / nrow(interaction.matrix)

data.plot1 <- data.frame(Interaction = "Pro Tumor", Frequency = interaction.freq.worse)
data.plot2 <- data.frame(Interaction = "Anti Tumor", Frequency = interaction.freq.better)
data.plot <- rbind(data.plot1, data.plot2)

box.plot <- ggplot(data = data.plot, aes(x = Interaction, y = Frequency)) + geom_boxplot() + theme_bw() + 
				geom_abline(slope = 0, intercept = expectation, col = "red", linetype = "dotted")

svg(file = output.file, height = 3, width = 2.5)
box.plot
dev.off()


#########obsolete###########z



#mutation.matrix <- clonality::create.mutation.matrix(maf, rem = FALSE)
#mutation.matrix <- SMDIC::maf2matrix(maf.input, percent = 0, nonsynonymous = TRUE)
#colnames(mutation.matrix) <- gsub("\\.", "-", colnames(mutation.matrix))
#mutation.matrix <- data.frame(mutation.matrix, check.names = F)



#variables <- cbind(binValues, bins)
#variables$samples <- gsub("\\w.\\w+.\\w+.\\w+$", "", variables$barcodes) %>% gsub("\\.", "-", .)

#gene.test.result <- list()
#for (gene in rownames(mutation.data)) {
#	data.test <- mutation.data[gene, ] %>% unlist %>% data.frame(samples = names(.), mutation = .)
#	data.test$mutation.stat <- ifelse(data.test$mutation == 0, 0, 1) ###binarize
#	data.test <- merge(variables, data.test, by.x = "samples")
#	
#	test.result <- glm(mutation.stat ~ bin1 + bin2 + bins, data = data.test, family = "binomial")
#	gene.test.result[[gene]] <- summary(test.result)
#	
#}