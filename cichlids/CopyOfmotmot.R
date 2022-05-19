#Load libraries
library(motmot)
library(tidyverse)
library(geiger)

#Load data for brain
path <- "cichlids/expression_data_and_trees/CichlidX_TPM_GeneExpressionMatrix_"
path2 <- "cichlids/expression_data_and_trees/"

br_data <- read.delim(paste0(path, "BR",".txt"))

#Annotate genes vs lnc
dictionary <- read.delim(paste0(path2, "GCF_001858045.1_ASM185804v2_genomic_gtf_gene.txt")) %>%
  select(geneID, biotype)

#Load tree
tree <- read.tree(paste0(path2, "intree")) %>% drop.tip("Calple")

#function to separate data into lnc or protein coding
separate_df <- function(matrix, type){
  step1 <- dictionary %>% filter(biotype == type)
  res <- matrix[which(row.names(matrix) %in% step1$geneID),]
  res
}

br_prot <- br_data %>% separate_df("protein_coding") %>% t() %>% treedata(tree, ., sort = TRUE)

#Now split the brain data into pieces
br_df <- br_prot$data 
br_motmot <- br_df[,5] %>% as.matrix()

#Find rate shifts
tm2 <- transformPhylo.ML(br_motmot, tree, model = "tm2", minCladeSize = 5, nSplits = 3)
sum_tm2 <- summary(tm2)
plot(sum_tm2)
cutOff <- calcCutOff(tree, n = 1000, mc.cores = 16, model = "tm2")
res2 <- summary(tm2, cutoff = cutOff)
saveRDS(res2, "cichlids/motmot2.rds")