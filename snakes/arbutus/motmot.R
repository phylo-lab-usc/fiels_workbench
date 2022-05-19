#Load libraries
library(motmot)
library(tidyverse)
library(geiger)

#Load data
tree <- read.nexus("snakes/data/Final_tree.nex")
data <- read.csv("snakes/data/vPhen_data.csv") %>% column_to_rownames(var = "species") %>%
  select(-c(Total, cm, Family, Author, Tech, Ref)) %>% as.matrix()

#Now split the brain data into pieces
motmot_df <- data[,"SVMP"] %>% as.matrix()

#Find rate shifts
tm2 <- transformPhylo.ML(motmot_df, tree, model = "tm2", minCladeSize = 3, nSplits = 3)
sum_tm2 <- summary(tm2)
plot(sum_tm2)
cutOff <- calcCutOff(tree, n = 1000, mc.cores = 4, model = "tm2")
res1 <- summary(tm2, cutoff = cutOff)
saveRDS(res1, "snakes/arbutus/motmot.rds")