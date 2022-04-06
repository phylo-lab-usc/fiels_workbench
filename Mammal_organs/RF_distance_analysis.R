#compare RF distance to phylosig

library(tidyverse)
library(phytools)
library(ape)
library(geiger)

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
  ggplot(aes(x = distance, y = log(phylosig))) + geom_point(position = position_jitter()) +
  theme_bw() + geom_smooth(method = "lm") + ylab("Log of Phylogenetic Half Life") + xlab("RF Distance") +
  ggtitle("Phylogenetic Inertia of Gene Trees as a function of RF Distance from Species Tree")

#Get phylogenetic K values
amniote_RPKM <- read.delim("Mammal_organs/Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.txt")

br_dat <- amniote_RPKM %>% select(hsa,contains(".br.")) 
br_avg_dat <- br_dat %>% group_by(hsa) %>% rename(Gene = hsa) %>% transmute("Homo sapiens" = mean(hsa.br.M.1, hsa.br.M.2, hsa.br.M.3, hsa.br.M.4, hsa.br.M.5, hsa.br.F.1),
                                                                            "Pan troglodytes" = mean(ptr.br.M.1, ptr.br.M.2, ptr.br.M.3, ptr.br.M.4, ptr.br.M.5, ptr.br.F.1),
                                                                            "Pan paniscus" = mean(ppa.br.M.1, ppa.br.F.2, ppa.br.F.1),
                                                                            "Gorilla gorilla" = mean(ggo.br.M.1, ggo.br.F.1),
                                                                            "Pongo pygmaeus" = mean(ppy.br.M.1, ppy.br.F.1),
                                                                            "Mus musculus" = mean(mmu.br.M.1, mmu.br.M.2, mmu.br.F.1),
                                                                            "Macaca mulatta" = mean(mml.br.F.1, mml.br.M.1, mml.br.M.2),
                                                                            "Monodelphis domestica" = mean(mdo.br.M.1, mdo.br.F.1),
                                                                            "Ornithorhynchus anatinus" = mean(oan.br.M.1, oan.br.F.1),
                                                                            "Gallus gallus" = mean(gga.br.M.1, gga.br.F.1))

br_avg_dat <- br_avg_dat %>% filter(Gene %in% names(gene_trees))
rm(br_dat)

format_expr_data <- function (avgdat) {
  temp <- avgdat %>% pull(Gene)
  avgdat <- avgdat %>% ungroup() %>% select(!Gene)
  dat <- flip(avgdat)
  colnames(dat) <- temp
  res <- dat %>% as.matrix()
  res
}

calculate_K_vals <- function(fmt_data, treelist){
  res <- vector(mode = "numeric", length = length(treelist))
  for(i in 1:length(treelist)){
    tmp <- tryCatch(phylosig(treelist[[i]], fmt_data[,i], test = TRUE), error = function(x)list(P = NA))
    res[i] <- tmp$P
  }
  res
}

Pvals <- br_avg_dat %>% format_expr_data() %>% calculate_K_vals(fmt_data = .,treelist = gene_trees)
RF_distance %>% mutate(Pvals) %>%
  ggplot(aes(x = distance, y = Pvals)) + geom_point(position = position_jitter()) +
  theme_bw() + geom_smooth(method = "lm") + ylab("Phylogenetic Signal P Value") + xlab("RF Distance") +
  ggtitle("Phylogenetic Signal of Gene Trees as a function of RF Distance from Species Tree")

RF_distance %>% mutate(Pvals) %>%
  ggplot(aes(group = distance, y = Pvals, fill = distance)) + geom_boxplot()
