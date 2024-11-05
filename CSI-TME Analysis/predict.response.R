library(dplyr); library(magrittr); 
library(glue); library(ggplot2)
#library(doParallel); library(doSNOW)
source('../functionsForInteractionAnalysis.R')

idh = "glioma"
cohort <- "immunotherapy"
dims = 10; method = "JADE"
remove.bcells = F

ggTheme <- theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1), #face = "bold"),
	axis.text.y = element_text(size = 12, face = "bold"),
	strip.background =element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 12),
	legend.title = element_blank(), legend.background = element_blank(), legend.position = "none",
	axis.title = element_text(size = 12), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	
##### load the interaction file#######

####load and bin the ICA components#####

####load patient clinical data#######
load(glue("ICA-{method}-{dims}.cellTypes.{cohort}.IDH{idh}-filtered-noNorm-optimalRank.Rda"))
interactionFile <- glue("ICA-JADE-{dims}.FDR0.25.robustInteractions-IDHmut-tcga.cgga-technical.txt")
interactions <- read.table(interactionFile, sep = "\t", header = T)
clinical <- readxl::read_excel("immunotherapy.glioma.clinical.xlsx", sheet = 1)

binMapList <- componentBinningFunction(pcaList = ica.res.list, factorization = "ICA", sample.names = NULL)
names(binMapList) <- gsub("\\.", "", names(binMapList))
interactions.robust <- interactions %>% subset(., FDR.x < 0.2 & pvalue.y < 0.05 & Reproduciblity > 7 & prognosis.x == prognosis.y)
#interactions.robust <- interactions %>% subset(., FDR.x < 0.2 & Reproduciblity > 7 & prognosis.x == prognosis.y)

if (remove.bcells == T) {
	indexes <- grep("Bcell", interactions.robust$string, invert = T)
	interactions.robust <- interactions.robust[indexes, ]
}


#interactions.robust <- interactions.robust[interactions.robust$prognosis.x == "Worst", ]
response.list <- list()
sample.matrix.list <- list();
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
	
	variables <- cbind(binValues, bins)
	
	clinicalData <- merge(clinical, variables, by = "barcodes") %>% na.omit()
	response.background <- 	table(clinicalData$responder.status)
	response.bins <- table(clinicalData$responder.status, clinicalData$bins)
	
	if (ncol(response.bins) == 2) {
		fraction.response <- response.bins['responder', '1'] / response.background['responder']
		fraction.noResponse <- response.bins['non-responder', '1'] / response.background['non-responder']
		
		################fisher p value##################
		fisher <- response.bins %>% .[ ,rev(colnames(.))] %>% fisher.test
		#if (prognosis  == "worst") {
		#	fisher <- response.bins %>% .[ ,rev(colnames(.))] %>% fisher.test ####columns are reversed to check enrichment of non-responders in 1 bin, i.e presence of interactions
		#}
		#if (prognosis  == "better") {
		#	fisher <- response.bins %>% fisher.test ####enrichment of 
		#}
		
		odds <- fisher$estimate; pvalue <- fisher$p.value
		response <- data.frame(string = interaction.string$string, Responder = fraction.response, "Non-Responder" = fraction.noResponse, odds = odds, pvalue = pvalue)
	} else {
		response <- data.frame(string = interaction.string$string, Responder = NA, "Non-Responder" = NA, odds = NA, pvalue = NA)
	}
	sample.matrix.list[[index]] <- bins %>% set_names(names(comp1.bins))
	
	response.list[[index]] <- response
}
response.df <- bind_rows(response.list)
max.finite <- sort(response.df$odds, decreasing  = T) %>% unique %>% .[2]
response.df$odds <- ifelse(is.finite(response.df$odds), response.df$odds, max.finite)

data.plot <- interactions.robust[ ,1:9]
colnames(data.plot) <- gsub("\\.x", "", colnames(data.plot))
data.plot <- merge(data.plot, response.df, by = "string")

#write.table(data.plot, glue("ICA-JADE-{dims}.immunotherapyResponse-technical.txt"), sep = "\t", quote = F, row.names = F)

##############plots################
box.plot <- ggplot(data = data.plot, aes(x = prognosis, y = odds)) + geom_boxplot() + theme_bw() + ggTheme + ggpubr::stat_compare_means() 
svglite::svglite(glue("plots/boxplot.immunotherapy.response.oddsRatio.svg"), height = 3, width = 3)
box.plot
dev.off()


data.plot <- reshape2::melt(data.plot, measure.vars = c("Responder", "Non.Responder"))
data.plot$interaction <- ifelse(data.plot$prognosis == "Better", "Anti-tumor", "Pro-tumor")
box.plot <- ggplot(data = data.plot, aes(x = variable, y = value)) + geom_boxplot() + theme_bw() + ggTheme + 
				facet_wrap(~interaction) + ggpubr::stat_compare_means() 


svglite::svglite(glue("plots/boxplot.immunotherapy.response.interactions.svg"), height = 3, width = 5)
box.plot
dev.off()

ggTheme <- theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1), #face = "bold"),
	axis.text.y = element_text(size = 12, face = "bold"),
	strip.background =element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 12),
	legend.title = element_blank(), legend.background = element_blank(), legend.position = "none",
	axis.title = element_text(size = 12), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	

################################################################################################################
################################################################################################################
################################################################################################################

col.annotations <- data.frame(response = clinical$responder.status) %>% set_rownames(clinical$barcodes)
row.annotations <- data.frame(prognosis = interactions.robust$prognosis.x) %>% set_rownames(interactions.robust$string)
sample.matrix <- do.call("rbind", sample.matrix.list) %>% set_rownames(interactions.robust$string)

svglite::svglite(glue("plots/heatmap.immunotherapy.response.interactions.svg"), height = 5, width = 3)
pheatmap::pheatmap(sample.matrix, annotation_col = col.annotations, annotation_row = row.annotations, clustering_method = "ward.D2", distance = "binary", show_rownames = F, show_colnames = F)
dev.off()



##########interaction load###########
sample.matrix <- do.call("rbind", sample.matrix.list) %>% set_rownames(interactions.robust$string)

int.pro <- interactions.robust[interactions.robust$prognosis.x == "Worst", 'string']
int.anti <- interactions.robust[interactions.robust$prognosis.x == "Better", 'string']
sample.matrix.pro <- sample.matrix[rownames(sample.matrix) %in% int.pro, ]
sample.matrix.anti <- sample.matrix[rownames(sample.matrix) %in% int.anti, ]


interaction.load <- colSums(sample.matrix.pro)
clinical.load <- merge(clinical, interaction.load, by.x = "barcodes", by.y = 0) %>% na.omit()
ggplot(data = clinical.load, aes(x = responder.status, y = y)) + geom_boxplot() + ggpubr::stat_compare_means()


interaction.load <- colSums(sample.matrix.anti)
clinical.load <- merge(clinical, interaction.load, by.x = "barcodes", by.y = 0) %>% na.omit()
ggplot(data = clinical.load, aes(x = responder.status, y = y)) + geom_boxplot() + ggpubr::stat_compare_means()




##################################Machine learning##################################
data.classification <- do.call("rbind", sample.matrix.list) %>% t() %>% data.frame()
index <- which(colSums(data.classification) == 0)
data.classification <- data.classification[ ,-index]
col.annotations <- na.omit(col.annotations)
data.classification <- merge(data.classification, col.annotations, by = 0) %>% set_rownames(.[ ,1]) %>% .[ ,-1]
#data.classification$response <- ifelse(data.classification$response == "non-responder", 1, 0)
data.classification$response <- factor(data.classification$response)


library(randomForest)
set.seed(123)

total.n <- nrow(data.classification)
train.prop <- 0.80
n <- floor(total.n * train.prop)

train.index <- sample(1:total.n, n, replace = F)
test.index <- setdiff(1:total.n, train.index)

train <- data.classification[train.index, ]
test <- data.classification[test.index, ]

rf_model <- randomForest(response ~ ., data = train, ntree = 1000, importance = T)                 
rf_model

predict_reg <- predict(logistic_model, test, type = "response")

					
logistic_model <- glm(response ~ ., data = train,family = "binomial")               
predict_reg <- predict(logistic_model, test, type = "response")

	
predict_reg <- ifelse(predict_reg > 0.5, 1, 0)
missing_classerr <- mean(predict_reg != test$response)
print(paste('Accuracy =', 1 - missing_classerr))


					