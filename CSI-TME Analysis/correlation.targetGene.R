library(dplyr); library(magrittr); library(glue)
library(doParallel); library(tidyr)
source('../functionsForInteractionAnalysis.R')

idh <- "mut"
target.gene <- "IDH1"

expression.file <- glue("../raw_gene_exp_matrix/tcga_idh_{idh}_tpm.txt")

expression.data <- read.table(expression.file, sep = "\t", header = T, row.names = 1)

log.tpm <- log(expression.data + 0.001) %>% data.matrix()
log.tpm <- log.tpm %>% data.matrix %>% preprocessCore::normalize.quantiles(.) %>% 
			set_rownames(rownames(log.tpm)) %>% set_colnames(colnames(log.tpm))
			
			
target.vector <- expression.data[target.gene, ] %>% t()

cor.output <- psych::corr.test(target.vector, y = t(expression.data), use = "pairwise", method = "spearman", adjust = "fdr", ci = F)

correlation <- cor.output$r %>% t() %>% reshape2::melt() %>% data.frame() %>% set_names(c("gene", "target", "correlation"))
pvalues <- cor.output$p %>% t() %>% reshape2::melt() %>% data.frame() %>% set_names(c("gene", "target", "pvalue"))
correlation$pvalue <- pvalues$pvalue
correlation$FDR <- p.adjust(correlation$pvalue, method = "fdr")
correlation$direction <- ifelse(correlation$correlation > 0, "positive", "negative")
correlation <- na.omit(correlation)
significant.correlation <- correlation %>% subset(., FDR < 0.05 & direction == "positive", select = 'gene')

write.table(correlation, "idh1.correlation.glioma.IDHmutant", sep = "\t", quote = F, row.names = F)