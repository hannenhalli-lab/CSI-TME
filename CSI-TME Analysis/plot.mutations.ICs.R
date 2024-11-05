			#####################frequency of genes with association################
			

library(dplyr); library(magrittr); library(glue)
library(tidyr); library(ggplot2); library(tidyverse)
source('../functionsForInteractionAnalysis.R')



idh = "mut"
cohort <- "tcga"
dims = 10; method = "JADE"
input.type = "interactions"
