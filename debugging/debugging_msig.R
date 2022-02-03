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
