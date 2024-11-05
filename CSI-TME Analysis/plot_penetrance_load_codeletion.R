#library(clusterProfiler)
library(dplyr); library(magrittr);
library(reshape2); library(ggplot2);
library(glue); 

source('../functionsForInteractionAnalysis.R')

library(survival)
library(survminer)

create_km_plot <- function(clinical.data, cancer_type, groups.data, variable.name, make_groups) {
	
	if (cancer_type == "use_all") {
		clinical.cancer <- clinical.data
	} else {
		clinical.cancer <- clinical.data[clinical.data$type %in% cancer_type, ]
	}
	clinical_matched <- merge(clinical.cancer, groups.data, by = "samples")
	clinical_matched$group <- clinical_matched[ ,variable.name]

	# survival analysis for KM plot
	fit <- survfit(Surv(time, status) ~ group, data = clinical_matched)

	levels.names <- fit$strata %>% names %>% gsub(".+=", "", .)

	p <- ggsurvplot(fit, data = clinical_matched, pval = TRUE, 
	              #title = paste(variable.name),
	              xlab = "Time (days)", ylab = "Survival Probability",
	              legend.title = NULL,
	              legend.labs = levels.names)
  
  return(p)
  
}

create.colors <- function(min.bound = 0, max.bound = 10, middle.value = 1, partition.length = 10) {
	#min.bound --> minimum value in the data to color for
	#max.bound --> maximum value in the data to color for
	#middle.value --> middle value to center around the colors
	#partition.length <- how many intervals each values should be partininted into below and aboce the middle.value
	
	
	by.lower <- middle.value / partition.length
	by.upper <- max.bound / partition.length

	breaks.lower <- seq(min.bound, middle.value, by = by.lower)
	breaks.upper <- seq(middle.value, max.bound, by = by.upper)[-1]

	breaks <- c(breaks.lower, breaks.upper)
	paletteLength <- length(breaks) - 1
	color.pallete <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
	return(list(colors = color.pallete, breaks = breaks))
	
}


ggTheme <- theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5), #face = "bold"),
	axis.text.y = element_text(size = 12, face = "bold"),
	strip.background =element_rect(fill = "white", color = "violetred4"), strip.text = element_text(colour = 'black', size = 12),
	legend.position = c(0.8, 0.4),
	axis.title = element_text(size = 12), panel.border = element_rect(color = "grey", fill = NA, size = 1))
	


idh <- "mut"; dims = 10; cohort = "CGGA693"
penetrance <- read.table(glue("CSIN_penetrance_codeletion_{cohort}.txt"), sep = "\t", header = T)
intLoad <- read.table(glue("CSIN_interactionLoad_TMEclass_{cohort}.txt"), sep = "\t", header = T)


if (cohort == "TCGA") {
	clinical_file = "AllAvialableClinicalDataMod.txt"
	clinical <- read.table(clinical_file, sep = "\t", header = T)
	coldel <- read.table("tcga_coldeletion.txt", sep = "\t", header = T)
	clinical <- clinical[clinical$type %in% c("GBM" ,"LGG"), ]
	clinical <- merge(clinical, coldel, by = "samples")
}
if (cohort == "CGGA693") {
	clinical_file = "../All_metadata-CCGA.csv"
	clinical <- read.csv(clinical_file,  header = T) %>% set_names(c('samples' ,'project', 'time', 'status', 'IDH', 'grade', 'age', 'subtype', 'codeletion', 'class'))
	clinical <- clinical[clinical$project == cohort, ]
}



data_plot <- penetrance[ ,c('string', 'IDH.O','IDH.A', 'prognosis')] %>% reshape2::melt()
svglite::svglite(file = glue("plots/boxplot_interaction_penetrance_codeletion_{cohort}.svg"), width = 2.5, height = 3)
ggplot(data = data_plot, aes(x = variable, y = value)) + geom_boxplot() + theme_bw() + ggTheme + 
	ggpubr::stat_compare_means(label = "p.format") + ylab("Interaction Penetrance")
dev.off()

data_plot <- penetrance[penetrance$dominance != "None", ]
fisher.matrix <- table(data_plot$prognosis, data_plot$dominance)

svglite::svglite(file = glue("plots/heatmap_codeletion_dominance_{cohort}.svg"), width = 3, height = 1.5)
pheatmap::pheatmap(fisher.matrix, cluster_row = F, cluster_col = F, display_numbers = T)
dev.off()


data_plot <- penetrance
svglite::svglite(file = glue("plots/scatterplot_penetrance_codeletion_{cohort}.svg"), width = 3, height = 3)
ggplot(data = data_plot, aes(x = IDH.A, y = IDH.O, col = dominance)) + theme_bw() + ylab("Oligodendroma") + xlab("Astrocytoma") +
	geom_point() + coord_cartesian(xlim = c(0, 0.5), ylim = c(0,0.5)) + geom_abline(slope = 1, intercept = 0) + ggTheme
dev.off()


#svglite::svglite(file = "plots/boxplot_penetrance_codeletion.svg", width = 5, height = 3)
#ggplot(data = penetrance, aes(x = dominance, y = HR)) + geom_boxplot() + facet_wrap(~codeletion) + theme_bw() + ggTheme + 
#	ggpubr::stat_compare_means(label = "p.format") + ylab("Interaction Penetrance")
#dev.off()

#####plot interaction load#####
data_plot <- intLoad %>% na.omit()
svglite::svglite(file = glue("plots/scatterplot_interactionLoad_{cohort}.svg"), width = 3, height = 3)
ggplot(data = data_plot, aes(x = Pro.tumor, y = Anti.tumor, col = TME_class)) + theme_bw() + ylab("Anti-tumor") + xlab("Pro-tumor") +
	geom_point() + coord_cartesian(xlim = c(0, 0.6), ylim = c(0,0.6)) + geom_abline(slope = 1, intercept = 0) + ggTheme
dev.off()


########KM plots

groups_data <- intLoad[ ,c('samples', 'class')] %>% na.omit() %>% set_colnames(c("samples", "class_idh"))
surv_plot_overall <- create_km_plot(clinical.data = clinical, cancer_type = "use_all", groups.data = groups_data, variable.name = "class_idh")


groups_data <- intLoad[ ,c('samples', 'class', 'TME_class')]  %>% na.omit()
groups_data$group <- paste(groups_data$class, groups_data$TME_class, sep = "_")
groups_data <- groups_data[groups_data$TME_class != 'None', ]
surv_plot_load <- create_km_plot(clinical.data = clinical, cancer_type = "use_all", groups.data = groups_data, variable.name = "group")



svglite::svglite(file = glue("plots/kmplot_codeletion_{cohort}.svg"), width = 4, height = 4)
surv_plot_overall$plot
dev.off()


svglite::svglite(file = glue("plots/kmplot_codeletion+TMEclass_{cohort}.svg"), width = 4, height = 4)
surv_plot_load$plot + guides(colour = guide_legend(ncol = 2))
dev.off()





######heatmap###
data_plot <- intLoad
data_plot <- data_plot[order(data_plot$TME_class, data_plot$oddsRatio), ] %>% set_rownames(NULL)
data_plot <- data_plot[ ,c('samples', 'Pro.tumor', 'Anti.tumor')]  %>% tibble::column_to_rownames(., var = "samples") %>% na.omit()

row_annotation_df <- intLoad[ ,c('samples', 'class', 'TME_class')] %>% tibble::column_to_rownames(., var = "samples") %>% na.omit()
colnames(row_annotation_df) <- c('subtype', 'TME_class')

ann_colors = list(
    TME_class = c("Pro-tumor" = "blue", "None" = "brown", "Anti-tumor" = "pink"),
    subtype = c("IDH-A" = "red", "IDH-O" = "green")
)


svglite::svglite(file = glue("plots/heatmap_TMEclass_load_{cohort}.svg"), width = 3, height = 6)
pheatmap::pheatmap(data_plot, annotation_row = row_annotation_df, show_rownames = F, cluster_rows = F, cluster_cols = F, annotation_colors = ann_colors)
dev.off()










######osbolete, not used #####
intLoad <- data.plot %>% group_by(samples, codeletion, int.type) %>% summarize(load = sum(value)/length(value)) %>% na.omit()
intLoad$codeletion <- ifelse(intLoad$codeletion == "Codel", "Oligodendroma", "Astrocytoma")

intLoad_astro <- intLoad[intLoad$codeletion == "Astrocytoma", ]
intLoad_oligo <- intLoad[intLoad$codeletion == "Oligodendroma", ]

intLoad_data <- reshape2::dcast(intLoad, formula = samples ~ int.type, value.var = "load") %>% tibble::column_to_rownames(., var = "samples")
intLoad_astro_data <- reshape2::dcast(intLoad_astro, formula = samples ~ int.type, value.var = "load") %>% tibble::column_to_rownames(., var = "samples")
intLoad_oligo_data <- reshape2::dcast(intLoad_oligo, formula = samples ~ int.type, value.var = "load") %>% tibble::column_to_rownames(., var = "samples")
intLoad_data <- intLoad_data[rowSums(intLoad_data) > 0, ]

my.pallete <- create.colors(min.bound = 0, max.bound = 0.5, middle.value = 0.1, partition.length = 10)

load_heatmap <- pheatmap(intLoad_astro_data, color = my.pallete[['colors']], breaks = my.pallete[['breaks']], clustering_method = "complete")
cluster_data <- cutree(load_heatmap$tree_row, k = 2) %>% data.frame(cluster = .)


intLoad_data <- intLoad_data[rownames(intLoad_data) %in% rownames(cluster_data),]
pheatmap(intLoad_data, color = my.pallete[['colors']], breaks = my.pallete[['breaks']], show_rownames = F, clustering_method = "complete", annotation_row = cluster_data)


gg1 <- merge(intLoad_data, cluster_data, by = 0)
boxplot(`Anti-tumor` ~ cluster, data = gg1)







