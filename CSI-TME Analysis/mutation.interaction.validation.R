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

mutationFile <- "CGGA.WEseq_286.20200506.txt"
mutation.data <- read.table(mutationFile, sep = "\t", header = T, row.names = 1, na.strings = "")
mutation.data[!is.na(mutation.data)] <- 1
mutation.data[is.na(mutation.data)] <- 0

mutation.data <- t(mutation.data)

idh = "mut"
cohort <- "cgga693"
dims = 10; method = "JADE"
remove.bcells = F

association.file <- glue("ICA-JADE-{dims}.FDR0.20.robustInteractions-IDH{idh}-association.Mutations.txt")
association.data <- read.table(association.file, sep = "\t", header = T)


#test.interactions <- association.data[association.data$prognosis == "Better" & association.data$FDR < 0.20, ]
test.interactions <- association.data[association.data$FDR < 0.20, ]
test.interactions <- test.interactions$predictor %>% unique()

#test.genes <- association.data[association.data$prognosis == "Better" & association.data$FDR < 0.20, ]
test.genes <- association.data[association.data$FDR < 0.20, ]
test.genes <- test.genes$response %>% unique()


####load patient clinical data#######
load(glue("ICA-{method}-{dims}.cellTypes.{cohort}.IDH{idh}-filtered-noNorm-optimalRank.Rda"))


binMapList <- componentBinningFunction(pcaList = ica.res.list, factorization = "ICA", sample.names = NULL)
names(binMapList) <- gsub("\\.", "", names(binMapList))



#interactions.robust <- interactions.robust[interactions.robust$prognosis.x == "Worst", ]
interaction.matrix.list <- list()
for (index in 1:length(test.interactions)) {
	interaction.string <- test.interactions[index]
	
	interaction.string <- strsplit(interaction.string, "_") %>% unlist()
	
	
	cellTypes <- strsplit(interaction.string[1], ":") %>% unlist()
	components <- strsplit(interaction.string[2], ":") %>% unlist()
	
	int.bin <- interaction.string[3]
	
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
colnames(interaction.matrix) <- test.interactions

#####filter mutation data#####


common.samples <- intersect(rownames(interaction.matrix), rownames(mutation.data))


mutation.data.filtered <- mutation.data[common.samples, test.genes]
mutation.freqeuncy <- apply(mutation.data, 1, function(x) sum(x > 0) / length(x))


interaction.matrix <- interaction.matrix[common.samples, ]


output.regression <- association.function(data1 = mutation.data.filtered, data2 = interaction.matrix, method.regression = "fisher")

output.regression <- bind_rows(output.regression, .id = "comparison")
output.regression$FDR <- p.adjust(output.regression$p.value, method = "fdr")


assocaition.validation <- output.regression[ ,c('comparison', 'Estimate', 'FDR')]



assocaition.comparison <- merge(association.data, assocaition.validation, by = "comparison")


assocaition.significant <- assocaition.comparison[assocaition.comparison$FDR.x < 0.2, ]

data.plot <- assocaition.significant[ ,c('Estimate.x', 'Estimate.y')] %>% set_colnames(c("Estimate.TCGA", "Estimate.CGGA"))


###replace Infinite estimates with max values for the cohort
data.plot <- apply(data.plot, 2, function(x) {ifelse(is.finite(x), x, max(x[is.finite(x)]))})
data.plot <- log2(data.plot + 0.001)

###########################Mutation load###############################
library(ggplot2)
output.file <- glue("plots/comparison.mutation.ionteraction.association.TCGA.CGGA.svg")

scatter.plot <- ggplot(data = data.plot, aes(x = Estimate.TCGA, y = Estimate.CGGA)) + geom_point() + geom_smooth(method = "lm") + 
						ggpubr::stat_cor(method = "spearman", col = "red") + theme_bw()

ggsave(file = output.file, height = 4, width = 4)
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