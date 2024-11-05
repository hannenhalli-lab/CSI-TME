#library(clusterProfiler)
library(dplyr); library(magrittr);
library(reshape2); library(ggplot2);
library(glue); library(ggpubr)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
	data.cor <- cbind(x, y) %>% na.omit() %>% data.frame()
    Cor <- abs(cor(data.cor$x, data.cor$y)) # Remove abs function if desired
    txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
    if(missing(cex.cor)) {
        cex.cor <- 0.4 / strwidth(txt)
    }
    text(0.5, 0.5, txt,
         cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}



#files <- list.files(path = ".", pattern = "survivalAnalysis.Genes-IDHmut") 
#files <- list.files(path = ".", pattern = "Malignant.txt") %>% grep("IDHwt", ., value = T)
#files <- list.files(path = ".", pattern = "survivalAnalysis.ICs-JADE-25") %>% grep(sampleType, ., value = T)
sampleType <- "mut"; dims = 10

files <- list.files(path = ".", pattern = "survivalAnalysis.ICs-JADE-10") %>% grep(sampleType, ., value = T)
file.names <- gsub(".+-|.txt", "", files)
names(files) <- file.names

survival.data <- lapply(files, read.table, sep = "\t", header = T)

survival.cross <- survival.data[[2]]
survival.same <- survival.data[c(1,3)] %>% purrr::reduce(left_join, by = c("cellType", "Component"))

svglite::svglite(file = glue("plots/HazardRatio.ICs-Cohorts-{sampleType}.svg"), width = 5, height = 5)
ggplot(data = survival.same, aes(x = HR.x, y =  HR.y)) + geom_point() + theme_bw() + geom_smooth(method = "lm") + facet_wrap(~cellType) + stat_cor(method = "spearman")
dev.off()


cellTypes <- survival.data[[1]]$cellType %>% unique()

data.cor.list <- list()
for (cellType in cellTypes) {
	discovery <- survival.data[['TCGA']] %>% .[grepl(cellType, .$cellType), ]
	validation <- survival.data[['cross']] %>% .[grepl(cellType, .$cellType), ]
	
	data.cor <- merge(discovery, validation, by = "Component")
	data.cor <- data.cor %>% group_by(cellType.y) %>% summarize(cor(HR.x, HR.y, method = "spearman")) %>% data.frame() %>% set_names(c("cellType", "correlation"))
	data.cor.list[[cellType]] <- data.cor
}

correlation.cross <- bind_rows(data.cor.list, .id = "cellType") %>% mutate(comparison = "across cell type")
correlation.same <- survival.same %>% group_by(cellType) %>% summarize(cor(HR.x, HR.y, method = "spearman")) %>% 
			data.frame() %>% set_names(c("cellType", "correlation")) %>% mutate(comparison = "same cell type")

data.plot <- rbind(correlation.same, correlation.cross)
write.table(data.plot, file = "correlation.survival.IC.cellTypes.txt", sep = "\t", quote = F, row.names = F)

data.plot$comparison <- factor(data.plot$comparison, levels = unique(data.plot$comparison))
svglite::svglite(file = glue("plots/HazardRatio.correlation.ICs-acrossCohorts-{sampleType}.svg"), width = 3.5, height = 3.5)
ggplot(data = data.plot, aes(x = comparison, y =  correlation)) + geom_boxplot() + theme_bw() +  stat_compare_means()
dev.off()
#########################
survival.data.HR <- survival.data[ ,c('HR.x', 'HR.y')] %>% set_names(file.names) 	  
svglite::svglite(file = glue("plots/HazardRatio.Genes-Cohorts-{sampleType}.svg"), width = 5, height = 5)
pairs(survival.data.HR, upper.panel = panel.cor,  lower.panel = panel.smooth)
dev.off()



##########For genes##########

files <- list.files(path = ".", pattern = "survivalAnalysis.Genes-IDHmut.*cellTypes.txt") 

survival.data <- lapply(files, read.table, sep = "\t", header = T)
names(survival.data) <- c("CGGA", "TCGA")
cellTypes <- survival.data[[1]]$cellType %>% unique()
cellTypeCombinations <- combn(cellTypes, 2) %>% t() %>% data.frame()
self <- data.frame(X1 = cellTypes, X2 = cellTypes)
cellTypeCombinations <- rbind(self, cellTypeCombinations)

data.cor.list <- list()
for (index in 1:nrow(cellTypeCombinations)) {
	cellType1 <- cellTypeCombinations[index, 'X1']
	cellType2 <- cellTypeCombinations[index, 'X2']
	
	discovery <- survival.data[['TCGA']] %>% .[grepl(cellType1, .$cellType), ]
	validation <- survival.data[['CGGA']] %>% .[grepl(cellType2, .$cellType), ]
	
	data.cor <- merge(discovery, validation, by = "Component")
	#data.cor <- data.cor %>% group_by(cellType.y) %>% summarize(cor(HR.x, HR.y, method = "spearman")) %>% data.frame() %>% set_names(c("cellType", "correlation"))
	correlation <- cor(data.cor$HR.x, data.cor$HR.y)
	cellType1 <- gsub("\\.txt", "", cellType1)
	cellType2 <- gsub("\\.txt", "", cellType2)
	
	data.cor.list[[index]] <- data.frame(cellType1 = cellType1, cellType2 = cellType2, correlation)
	
}


data.correlation <- do.call("rbind", data.cor.list)
data.correlation$type <- ifelse(data.correlation$cellType1 ==  data.correlation$cellType2, "same cell type", "across cell types")
write.table(data.correlation, file = "correlation.survival.genes.cellTypes.txt", sep = "\t", quote = F, row.names = F)


