#compare RF distance to phylosig

species_tree <- read.tree("Mammal_organs/species_phylogeny/species_names.nwk") %>%
  drop.tip("Pan_paniscus")
species_tree$tip.label <- species_tree$tip.label %>%
  gsub(pattern = "_", replacement = " ", x = .)
gene_trees <- readRDS("Mammal_organs/gene_family_phylogenies/genefamilytrees")

tip.compare <- function(treelist){
  res <- vector(mode = "numeric", length = length(treelist))
  for(i in 1:length(treelist)){
    genetr <- treelist[[i]]
    tmp <- keep.tip(species_tree, genetr$tip.label)
    res[[i]] <- dist.topo(genetr, tmp)
  }
  res
}

RF_distance <- tip.compare(gene_trees) %>% as.data.frame() %>% rename(distance = ".")

RF_distance %>%
  ggplot(aes(x = distance)) + geom_bar() + theme_bw() + scale_x_binned()

get_alpha <- function(fit_obj){
  alpha <- fit_obj$opt$alpha
  res <- log(2)/alpha
  res
}

OUfits <- readRDS("Mammal_organs/gene_family_phylogenies/arbutus/justOU/fit_br")
alpha_list <- map_dbl(OUfits, get_alpha)

RF_distance <- RF_distance %>% mutate(phylosig = alpha_list)

RF_distance %>%
  ggplot(aes(x = distance, y = phylosig)) + geom_point() + theme_bw()
