#debugging m.sig

library(geiger)
library(ape)
library(arbutus)
library(OUwie)
library(tidyverse)

#Create and modify phylogenetic tree
simtree <- sim.bdtree()
simtree$node.label <- rep(1,99)
dat <- data.frame(simtree$tip.label) %>% mutate(Reg = 1) %>% rename(Genus_species = simtree.tip.label)

#Simulate data
simulate_data <- function ( phy, df ){
  sim <- OUwie.sim(phy, df, alpha = c(1.0, 1.0), sigma.sq = c(0.9, 0.9), theta0 = 1, theta = c(1.0, 1.0))
  rownames(sim) <- sim$Genus_species
  fitdf <- sim %>% select(X)
  fit <- fitContinuous(phy, fitdf, model = "OU")
  list(fit)
}
sim_fivehundred <- replicate(n = 500, expr = simulate_data(simtree, dat))

#Step by step, run arbutus as a whole
unit_trees <- map(sim_fivehundred, make_unit_tree)
obs <- calculate_pic_stat(unit_trees)
sim.dat <- simulate_char_unit(unit_trees)
sim <- calculate_pic_stat(sim.dat)
res <- compare_pic_stat(obs, sim)
arbutus::pic_stat_msig_plot(unit_trees[[5]])
plot(res)

#Now, run arbutus on each fit individually
arbutus_fivhund <- sim_fivehundred %>% map(arbutus)
pvals <- arbutus_fivhund %>% map_df(pvalue_arbutus)
pvals %>% ggplot(aes(x = m.sig)) + geom_histogram() + theme_bw() + ggtitle("M.sig p values") + xlab("P-values")

#Try running with different trees each time?
whole_process <- function(){
  tr <- sim.bdtree()
  tr$node.label <- rep(1,99)
  df <- data.frame(tr$tip.label) %>% mutate(Reg = 1) %>% rename(Genus_species = tr.tip.label)
  fitted_sim <- simulate_data(tr, df) 
  res <- fitted_sim %>% map(arbutus)
  res
}

dif_trees <- replicate(100, whole_process())
pvals_dt <- dif_trees %>% map_df(pvalue_arbutus)
pvals_dt %>% ggplot(aes(x = m.sig)) + geom_histogram() + theme_bw() + ggtitle("M.sig p values") + xlab("P-values")
piv_dt <- pvals_dt %>% pivot_longer(cols = everything(), values_to = "pvals")
piv_dt %>% ggplot(aes(x = pvals)) + geom_histogram() + theme_bw() + facet_wrap(~name) +
  ggtitle("p values") + xlab("P-values")

#Try with inadequate model
simulate_wrong <- function( phy, df ){
  sim <- OUwie.sim(phy, df, alpha = c(1.0, 1.0), sigma.sq = c(0.9, 0.9), theta0 = 1, theta = c(1.0, 1.0))
  rownames(sim) <- sim$Genus_species
  fitdf <- sim %>% select(X)
  fit <- fitContinuous(phy, fitdf, model = "EB")
  list(fit)
}

wrong_process <- function(){
  tr <- sim.bdtree()
  tr$node.label <- rep(1,99)
  df <- data.frame(tr$tip.label) %>% mutate(Reg = 1) %>% rename(Genus_species = tr.tip.label)
  fitted_sim <- replicate(500, simulate_wrong(tr, df)) 
}

step_by_step <- function(fitlist){
  unit_trees <- map(fitlist, make_unit_tree)
  obs <- calculate_pic_stat(unit_trees)
  sim.dat <- simulate_char_unit(unit_trees)
  sim <- calculate_pic_stat(sim.dat)
  res <- compare_pic_stat(obs, sim)
  res
}

inad_fit <- wrong_process()
inad_arb <- step_by_step(inad_fit)
plot(inad_arb)

## Get plots with real data
#Load newick format phylogeny and convert to something recognizable by ape/geiger
species_phylo <- read.tree(file = "Mammal_organs/species_phylogeny/species_names.nwk")

#Rename tips
species_phylo$tip.label <- gsub("_", " ", species_phylo$tip.label)

amniote_RPKM <- read.delim("Mammal_organs/Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.txt")

br_avg_dat <- amniote_RPKM %>% select(hsa,contains(".br.")) %>% group_by(hsa) %>% rename(Gene = hsa) %>% transmute("Homo sapiens" = mean(hsa.br.M.1, hsa.br.M.2, hsa.br.M.3, hsa.br.M.4, hsa.br.M.5, hsa.br.F.1),
                                                                            "Pan troglodytes" = mean(ptr.br.M.1, ptr.br.M.2, ptr.br.M.3, ptr.br.M.4, ptr.br.M.5, ptr.br.F.1),
                                                                            "Pan paniscus" = mean(ppa.br.M.1, ppa.br.F.2, ppa.br.F.1),
                                                                            "Gorilla gorilla" = mean(ggo.br.M.1, ggo.br.F.1),
                                                                            "Pongo pygmaeus" = mean(ppy.br.M.1, ppy.br.F.1),
                                                                            "Mus musculus" = mean(mmu.br.M.1, mmu.br.M.2, mmu.br.F.1),
                                                                            "Macaca mulatta" = mean(mml.br.F.1, mml.br.M.1, mml.br.M.2),
                                                                            "Monodelphis domestica" = mean(mdo.br.M.1, mdo.br.F.1),
                                                                            "Ornithorhynchus anatinus" = mean(oan.br.M.1, oan.br.F.1),
                                                                            "Gallus gallus" = mean(gga.br.M.1, gga.br.F.1)) %>%
  ungroup() %>% slice_head(n = 1000)

format_expr_data <- function (avgdat) {
  temp <- avgdat %>% pull(Gene)
  avgdat <- avgdat %>% ungroup() %>% select(!Gene)
  dat <- flip(avgdat)
  colnames(dat) <- temp
  dat
}

#Running fitcontinuous

runFC <- function (dat){
  fitResults <- vector(mode = "list", length = ncol(dat))
  tdf <- treedata(species_phylo, dat, sort = TRUE)
  phy <- tdf$phy
  data <- tdf$data
  for(j in 1:ncol(dat)){
    fitBM <- fitContinuous(phy, data[,j], model = "BM")
    fitOU <- fitContinuous(phy, data[,j], model = "OU")
    fitEB <- fitContinuous(phy, data[,j], model = "EB")
    aic <- c(fitBM$opt[["aic"]], fitOU$opt[["aic"]], fitEB$opt[["aic"]])
    fit <- ifelse(min(aic) == aic[1], list(c(fitBM, model = "BM")), 
                  ifelse(min(aic) == aic[2], list(c(fitOU, model = "OU")), 
                         list(c(fitEB, model = "EB"))))
    fitResults[j] <- fit
  }
  fitResults
}

total_process <- function (avgdat, part){
  exp <- format_expr_data(avgdat)
  fit <- runFC(exp)
  fit_name <- paste0("debugging/fit_", part)
  saveRDS(fit, file = fit_name)
  arb <- step_by_step(fit)
  arb_name <- paste0("debugging/arb_", part)
  saveRDS(arb, file = arb_name)
  plot(arb_name)
  ggsave("debugging/br_hist.png")
}


