
library(dplyr); library(magrittr);
library(glue); library(ggplot2)
library(ggpubr)

source('../../melanomaSplicing/generalFunctions.R')


idh = "mut";
cohort <- "TCGA"
dims = "10"
geneset.file <- "geneSets.Compilation.xlsx"


load(glue("ICA-JADE-{dims}.cellTypes.{tolower(cohort)}.IDH{idh}-filtered-noNorm-optimalRank.Rda"))
load(glue("ICA-JADE-{dims}.signatureGenes.{tolower(cohort)}.IDH{idh}--filtered-noNorm-optimalRank.Rda"))



readxl::excel_sheets(geneset.file)
#sheets <- c("Venteicher"); cellType <- "Malignant"
#sheets <- c("misgdb"); cellType <- "Tcell"
#sheets <- c("Wechter.TimeCourse"); cellType <- "Tcell"
#sheets <- c("CellAge"); cellType <- "Tcell"
sheets <-c("Biomarkers"); cellType <- "Tcell"
#sheets <-c("SnG.Clusters"); cellType <- "Tcell"
#sheets <- c("PSS"); cellType <- "Malignant"

#sheets <- readxl::excel_sheets(geneset.file)
#gene.sets <- lapply(sheets, function(x) readxl::read_excel(geneset.file, sheet = x))
#names(gene.sets) <- sheets

gene.set <- readxl::read_excel(geneset.file, sheet = sheets)
gene.set <- gene.set[ ,c('Gene', 'State')]
categories <- unique(gene.set$State)


ica.res <- ica.res.list[[cellType]]
S <- ica.res$S

S <- data.frame(Gene = rownames(S), S)
S.background <- S %>% mutate(State = "All Genes")
S.signatures <- merge(S, gene.set, by = "Gene")

S.combine <- rbind(S.signatures, S.background) %>% reshape2::melt()

#data.plot <- S.combine %>% subset(., State %in% c("Stemness program",  "All Genes") & variable == "IC.7")
#data.plot <- S.combine %>% subset(., State %in% c("Biomarker",  "All Genes") & variable == "IC.1")
data.plot <- S.combine %>% subset(., variable == "IC.7")

####for stemmness and CD44
output.file = glue("plots/boxplot.loadings.PSS.IC7.{cellType}.svg")
comparisons <- list(c("All Genes", "PSS1"), c("All Genes", "PSS2"))
svglite::svglite(file = output.file, height = 3, width = 3)
ggplot(data = data.plot, aes(x  = State, y = value)) + geom_boxplot() + stat_compare_means(label = "p.format", comparisons = comparisons) + 
	geom_point(data = data.plot[data.plot$Gene == "CD44", ], aes(x = "PSS1", y = value, col = "red")) +
	theme_bw() + ylab("Gene Weights IC7") + facet_wrap(~variable)
dev.off()




#box.plot <- ggplot(data = data.plot, aes(x  = variable, y = value)) + geom_boxplot() + stat_compare_means(label = "p.format", label.x = 1.5) + 
#	theme_bw() + ylab("Gene Weights IC1") + facet_wrap(~State)

explore = T
	
if (explore) {
	data.plot <- S.combine 
	box.plot <- ggplot(data = data.plot, aes(x  = State, y = value)) + geom_boxplot() + facet_wrap(~variable) +
	stat_compare_means(label = "p.format", label.x = 1.5, label.y = 3) + theme_bw() + ylab("Gene Weights IC1")
	
	svglite::svglite(file = glue("plots/boxplot.loadings.{sheets}.allICs.{cellType}.svg"), height = 5, width = 4.5)
	box.plot
	dev.off()
		
}

svglite::svglite(file = glue("plots/boxplot.loadings.stemness.IC7.{cellType}.svg"), height = 3, width = 3)
box.plot
dev.off()



###########select markers##########

cellType <- "Malignant"
ica.res <- ica.res.list[[cellType]]
S <- ica.res$S
S <- data.frame(Gene = rownames(S), S)

select.marker <- "CD44"
select.expression <- S[select.marker, 2:11] %>% unlist()
order.index <- paste('IC', 1:10, sep = ".")
select.expression <- select.expression[order.index]
select.expression <- data.frame(IC = names(select.expression), Weights = select.expression)

bar.plot <- ggplot(data = select.expression, aes(x = IC, y = Weights, group = 1)) + geom_bar(stat = "identity", fill = "white", col = "black") + theme_bw()
svglite::svglite(file = glue("plots/barplot.loadings.CD44.ICs.{cellType}.svg"), height = 3, width = 5)
bar.plot
dev.off()


####################################################
####################################################

cellType <- "Tcell"
ica.res <- ica.res.list[[cellType]]
S <- ica.res$S
S <- data.frame(Gene = rownames(S), S)

#select.marker <- c("CD57", "KLRG1") #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5578132/
#select.marker <- c("CD27", "CD28") #####loss of these markers in sensence 
#select.marker <- c("SESN3", "FOXO3", "STK17A")
#select.marker <- "CDKN1A"
select.marker <- "CD4"
select.expression <- S[rownames(S) %in% select.marker, 2:11] %>% unlist()
order.index <- paste('IC', 1:10, sep = ".")
select.expression <- select.expression[order.index]
select.expression <- data.frame(IC = names(select.expression), Weights = select.expression)

bar.plot <- ggplot(data = select.expression, aes(x = IC, y = Weights, group = 1)) + geom_bar(stat = "identity", fill = "white", col = "black") + theme_bw()
#svglite::svglite(file = glue("plots/barplot.loadings.CD44.ICs.{cellType}.svg"), height = 3, width = 5)
bar.plot
dev.off()








######obsolete#######
par(mfrow = c(2, 5))
for (index in 1:10) {
S.s <- S[ ,index]
gg1 <- S.s[names(S.s) %in% set.genes]
boxplot(gg1, S.s)

}