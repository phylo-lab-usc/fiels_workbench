#Parse Gene Trees into list

library(ape)
library(parallel)

tree_names <- list.files("Mammal_organs/gene_family_phylogenies/trees/") %>% as.list()

parse_trees <- function ( treestring ){
  tree <- read.tree(file = paste0("Mammal_organs/gene_family_phylogenies/trees/", treestring))
  tree
}

treelist <- mclapply(tree_names, parse_trees, mc.cores = 12)

saveRDS(treelist, file = "Mammal_organs/gene_family_phylogenies/genetrees")