#Multirate arbutus test 3: Test how regime affects adequacy

library(tidyverse)
library(ape)
library(geiger)
library(OUwie)
library(broom)
library(arbutus)

#Now, simulate birth-death tree

tr <- sim.bdtree(n = 128)

dat1 <- data.frame(tr$tip.label) %>% mutate(Reg = c(rep(1, 64), rep(2, 64))) %>% rename(Genus_species = tr.tip.label)
dat2 <- data.frame(tr$tip.label) %>% mutate(Reg = c(rep(1, 64), rep(2, 34), rep(3,30))) %>% rename(Genus_species = tr.tip.label)

#First test what happens when fitted model has more and less regimes than true model

#Use different trees for regular OU and OU variants
tree_OUwie <- tr
tree_OUwie$node.label <- c(rep(1, 64), rep(2, 63))

tree_false <- tr
tree_false$node.label <- c(rep(1, 64), rep(2, 33), rep(3,30))

sim_and_fit_arbutus <- function (tree, dat) {
  df <- OUwie.sim(tree, dat, alpha = c(1e-10, 1e-10), sigma.sq = c(0.45, 0.90), theta0 = 1.0, theta = c(0, 0))
  df_fix <- df
  row.names(df_fix) <- df_fix$Genus_species
  df_fix <- df_fix %>% select(X)
  df_false <- dat2
  df_false$X <- df$X
  bms <- OUwie(tree, df, model = "BMS") %>% arbutus() %>% pvalue_arbutus() #True model, 2 regimes
  bms2 <- OUwie(tree_false, df_false, model = "BMS") %>% arbutus() %>% pvalue_arbutus() #False model, 3 regimes
  bm <- fitContinuous(tree, df_fix, model = "BM") %>% arbutus() %>% pvalue_arbutus() #BM, 1 regime
  res <- list(bms = bms, bm = bm, bms2 = bms2)
  res
}

#run the sims
adequacy <- replicate(1000, sim_and_fit_arbutus(tree_OUwie, dat1))

#retrieve pvals
bms <- adequacy["bms",] %>% map_df(function(x)x) %>%
  pivot_longer(cols = everything(), names_to = "statistic", values_to = "pvalue") %>%
  mutate(model = "bms")

bm <- adequacy["bm",] %>% map_df(function(x)x) %>%
  pivot_longer(cols = everything(), names_to = "statistic", values_to = "pvalue") %>%
  mutate(model = "bm")

bms2 <- adequacy["bms2",] %>% map_df(function(x)x) %>%
  pivot_longer(cols = everything(), names_to = "statistic", values_to = "pvalue") %>%
  mutate(model = "bms_false")

#fuse
fuse_df <- full_join(bms, bm) %>% full_join(bms2) 
saveRDS(fuse_df,"Arbutus_Exploration/RDSfiles/multirate_data_v3.rds")

#plot
fuse_df %>% 
  ggplot(aes(y = pvalue, x = model, fill = (model))) + geom_violin() + geom_boxplot(width = 0.5) + facet_wrap(~statistic) + theme_bw()
ggsave("Arbutus_Exploration/Figures/multirate_regime_study.png")