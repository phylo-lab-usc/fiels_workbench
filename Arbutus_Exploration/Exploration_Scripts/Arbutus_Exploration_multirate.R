#Multirate arbutus test

library(tidyverse)
library(ape)
library(geiger)
library(OUwie)
library(broom)
library(arbutus)

#Now, simulate birth-death tree

tr <- sim.bdtree(n = 128)

dat1 <- data.frame(tr$tip.label) %>% mutate(Reg = c(rep(1, 64), rep(2, 64))) %>% rename(Genus_species = tr.tip.label)
dat2 <- data.frame(tr$tip.label) %>% mutate(Reg = 1) %>% rename(Genus_species = tr.tip.label)

#Use different trees for regular OU and OU variants
tree_OUwie <- tr
tree_OUwie$node.label <- c(rep(1, 64), rep(2, 63))

tree_OU <- tr
tree_OU$node.label <- rep(1,127)

sim_and_fit_arbutus <- function (tree, dat) {
  df <- OUwie.sim(tree, dat, alpha = c(1e-10, 1e-10), sigma.sq = c(0.45, 0.90), theta0 = 1.0, theta = c(0, 0))
  df_fix <- df
  row.names(df_fix) <- df_fix$Genus_species
  df_fix <- df_fix %>% select(X)
  bms <- OUwie(tree, df, model = "BMS") %>% arbutus() %>% pvalue_arbutus()
  bm <- fitContinuous(tree, df_fix, model = "BM") %>% arbutus() %>% pvalue_arbutus()
  ou <- fitContinuous(tree, df_fix, model = "OU") %>% arbutus() %>% pvalue_arbutus()
  eb <- fitContinuous(tree, df_fix, model = "EB") %>% arbutus() %>% pvalue_arbutus()
  res <- list(bms = bms, bm = bm, ou = ou, eb = eb)
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

ou <- adequacy["ou",] %>% map_df(function(x)x) %>%
  pivot_longer(cols = everything(), names_to = "statistic", values_to = "pvalue") %>%
  mutate(model = "ou")

eb <- adequacy["eb",] %>% map_df(function(x)x) %>%
  pivot_longer(cols = everything(), names_to = "statistic", values_to = "pvalue") %>%
  mutate(model = "eb")

#fuse
fuse_df <- full_join(bms, bm) %>% full_join(ou) %>% full_join(eb)
saveRDS(fuse_df,"Arbutus_Exploration/RDSfiles/multirate_data.rds")

#plot
fuse_df %>% 
  ggplot(aes(y = pvalue, x = model, fill = (model))) + geom_violin() + geom_boxplot(width = 0.5) + facet_wrap(~statistic) + theme_bw()
ggsave("Arbutus_Exploration/Figures/multirate_comparison")