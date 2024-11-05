#####in this script, we test the differential association between mutation and interactions between early vs. late samples

####inputs:
####clinical data with information on tumor grades######
####mutation matrix######
####interaction matrix#####
####significant mutation ~ interaction associations#####

library(dplyr)
library(glue)
library(magrittr)
library(simsalapar)

cancerType <- c("GBM", "LGG")
clinical <- read.table("/Users/singha30/Desktop/AllAvialableClinicalDataMods.txt", sep = "\t", header = T)
clinical <- clinical[clinical$type %in% cancerType, ]

codel <- read.csv("../All_metadata-CCGA.csv", header = T) %>% set_names(c('samples' ,'project', 'time', 'status', 'IDH', 'grade', 'age', 'subtype', 'codel'))
codel <- codel[grep("TCGA", codel$samples), c('IDH', 'samples', 'codel', 'grade')]
codel$samples <- gsub("\\.", "-", codel$samples) %>% gsub("-\\w+-\\w+-\\w+-\\w+$", "", .)
clinical <- merge(clinical, codel, by = "samples")

clinical <- clinical[clinical$IDH == "Mut", ]

grade.data <- clinical[ ,c('samples', 'grade')] %>% na.omit()
grade.data$grade.level <- ifelse(grade.data$grade == "G2", "Low", "High")
grade.data$grade <- factor(grade.data$grade)

mutation.matrix <- read.table("mutation.matrix.tcga.txt", sep = "\t", header = T, row.names = 1)
interaction.matrix <- read.table("interaction.matrix.tcga.txt", sep = "\t", header = T, row.names = 1)


mutation.interaction <- read.table("ICA-JADE-10.FDR0.20.robustInteractions-IDHmut-association.Mutations.txt", sep = "\t", header = T)
mutation.interaction <- mutation.interaction[mutation.interaction$FDR < 0.2, ]


rownames(mutation.matrix) <- gsub("-\\w+$", "", rownames(mutation.matrix))
rownames(interaction.matrix) <- gsub("-\\w+$", "", rownames(interaction.matrix))


mutation.matrix <- data.matrix(mutation.matrix)
interaction.matrix <- data.matrix(interaction.matrix)

variable.combniations <- unique(mutation.interaction$comparison) %>% unique()


estimate.association.grades.list <- list()
ratio.interaction.grades.list <- list()
for (variable.combniation in variable.combniations) {
	variable.partners <- strsplit(variable.combniation, "-") %>% unlist()
	gene <- variable.partners[1]
	interaction <- variable.partners[2]
	interaction <- gsub(":", ".", interaction) 
	
	samples.mutation <- mutation.matrix[ ,gene] %>% data.frame(samples = names(.), mutation = .) 
	samples.interaction <- interaction.matrix[ ,interaction] %>% data.frame(samples = names(.), interaction = .)
	
	combined.data <- merge(samples.mutation, samples.interaction, by = "samples")
	
	combined.data <- merge(combined.data, grade.data, by = "samples")
	combined.data$membership <- ifelse(combined.data$mutation == 1 & combined.data$interaction == 1, 1, 0)
	
	association <- tryCatch.W.E(glm(grade ~ membership + mutation + interaction, data = combined.data, family = "binomial"))
	#association <- tryCatch.W.E(glm(grade ~ membership, data = combined.data, family = "binomial"))
	coefficents.association <- summary(association$value)[['coefficients']] %>% data.frame
	colnames(coefficents.association) <- c("Estimate.Grade", "StdErr", "z.value", "p.value")
	warning <- association$warning
	
	if (is.null(warning)) {
		warning <- "No"
	} else {
		warning <- "Yes"
	}
	
	coefficents.association$warning <- warning
	estimate.association.grades.list[[variable.combniation]] <- coefficents.association
	
	combined.data.mutant <- combined.data[combined.data$mutation == 1, ]
	combined.data.wildtype <- combined.data[combined.data$mutation == 0, ]
	
	if (any(combined.data.mutant$interaction == 1) & any(combined.data.wildtype$interaction == 1)) {
		stat.mutant <- table(combined.data.mutant$interaction, combined.data.mutant$grade.level) %>% data.matrix
		stat.wildtype <- table(combined.data.wildtype$interaction, combined.data.wildtype$grade.level) %>% data.matrix
	
	
		stat.mutant <- stat.mutant['1', ]
		stat.wildtype <- stat.wildtype['1', ]
	
		stat.matrix <- matrix(nrow = 2, ncol = 2) %>% set_colnames(c("Low", "High")) %>% set_rownames(c("Mutant", "WildType"))
		stat.matrix[1, ] <- c(stat.mutant['Low'], stat.mutant['High'])
		stat.matrix[2, ] <- c(stat.wildtype['Low'], stat.wildtype['High'])
	
		ratio.mutant <- stat.matrix['Mutant', 'Low'] / stat.matrix['Mutant','High']
		ratio.wildtype <- stat.matrix['WildType', 'Low'] / stat.matrix['WildType','High']
	
	
		ratio.interaction.grades <- data.frame(mutant = ratio.mutant, wildtype = ratio.wildtype, odds = ratio.mutant / ratio.wildtype)
		ratio.interaction.grades.list[[variable.combniation]] <- ratio.interaction.grades
	}
	
	
	#association.mutant <- tryCatch.W.E(glm(grade ~ interaction, data = combined.data.mutant, family = "binomial"))
	#coefficents.mutant <- summary(association.mutant$value)[['coefficients']] %>% data.frame
	#colnames(coefficents.mutant) <- c("Estimate", "StdError", "t.value", "p.value")
	
	
	#association.wildtype <- tryCatch.W.E(glm(grade ~ interaction, data = combined.data.wildtype, family = "binomial"))
	#coefficents.wildtype <- summary(association.wildtype$value)[['coefficients']] %>% data.frame
	#colnames(coefficents.wildtype) <- c("Estimate", "StdError", "t.value", "p.value")
	
}

estimate.association.grades <- bind_rows(estimate.association.grades.list, .id = "comparison")
indexes <- grep("membership", rownames(estimate.association.grades))

estimate.association.grades <- estimate.association.grades[indexes, ]


estimate.association.grades <- merge(estimate.association.grades, mutation.interaction, by = "comparison")

write.table(estimate.association.grades, file = "TCGA.assocation.mutation.interactions.earlyVsLate.txt", sep = "\t", quote = F, row.names = F)


library(ggplot2); library(ggpubr)

ggplot(data = estimate.association.grades, aes(x = prognosis, y = Estimate.Grade)) + geom_boxplot() + stat_compare_means(label = "p.format", label.x = 1.5, label.y = 2) + coord_cartesian(ylim = c(-2.5, 2.5)) + theme_bw()

ggsave("plots/effect.tumorGrade.mutationInteractionAssociation.svg", height = 3, width = 2.5)
dev.off()



estimate.association.grades <- bind_rows(ratio.interaction.grades.list, .id = "comparison")
estimate.association.grades <- merge(estimate.association.grades, mutation.interaction, by = "comparison")

estimate.association.grades <- estimate.association.grades[ ,c('comparison', 'mutant', 'wildtype', 'odds', 'prognosis')]
estimate.association.grades <- reshape2::melt(estimate.association.grades) 
ggplot(estimate.association.grades, aes(x = variable, y = value)) + facet_wrap(~prognosis) + geom_boxplot() + 
	stat_compare_means(label = "p.format", label.x = 1.5, label.y = 2) + coord_cartesian(ylim = c(0, 2)) + theme_bw()
ggsave("plots/ratio.tumorGrade.mutationInteractionAssociation.svg", height = 3, width = 4)
dev.off()	
	
estimate.association.grades <- estimate.association.grades[ ,c('comparison', 'mutant', 'wildtype', 'odds', 'prognosis')]
ggplot(estimate.association.grades, aes(x = prognosis, y = odds)) + geom_boxplot() + geom_hline(yintercept = 1, lty = 2, col = "red") + ggpubr::stat_compare_means(label = "p.format", label.x = 1.5, label.y = 3.5) + coord_cartesian(ylim = c(0, 4)) + theme_bw()
ggsave("plots/oddsRatio.tumorGrade.mutationInteractionAssociation.svg", height = 3, width = 2.5)
dev.off()

	
	####scatter plot###
