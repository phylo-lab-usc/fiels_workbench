#Arbutus exploration continued

#First, load required packages

library(tidyverse)
library(ape)
library(geiger)
library(OUwie)
library(broom)
library(arbutus)


#Now, simulate birth-death tree

tr <- sim.bdtree(n = 128)

dat1 <- data.frame(tr$tip.label) %>% mutate(Reg = 1) %>% rename(Genus_species = tr.tip.label)
tree_OUwie <- tr
tree_OUwie$node.label <- rep(1, 127)

sim_and_fit_arbutus2 <- function (tree, dat, model) {
  ifelse(model == "OU", {
    df <- OUwie.sim(tree, data = dat, alpha = c(1.0, 1.0), sigma.sq = c(0.9, 0.9), theta0 = 1, theta = c(1.0, 1.0))
  }, 
  ifelse(model == "MV", {
    df <- OUwie.sim(tree, dat, alpha = c(1.0, 1.0), sigma.sq = c(0.45, 0.9), theta0 = 1.0, theta = c(1.0, 2.0))
  },
  ifelse(model == "MA", {
    df <- OUwie.sim(tree, dat, alpha = c(1.0, 0.5), sigma.sq = c(0.45, 0.45), theta0 = 1.0, theta = c(1.0, 2.0))
  },
  if(model == "MVA") df <- OUwie.sim(tree, dat, alpha = c(1.0, 0.5), sigma.sq = c(0.45, 0.9), theta0 = 1.0, theta = c(1.0, 2.0)))))
  
  df_fix <- df
  row.names(df_fix) <- df_fix$Genus_species
  df_fix <- df_fix %>% select(X)
  fit <- fitContinuous(tree, df_fix, model = "OU")
  a <- arbutus(fit)
  a
}

#run the sims
OU_adequacy <- replicate(1000, sim_and_fit_arbutus2(tree_OUwie, dat1, model = "OU"))
MV_adequacy <- replicate(1000, sim_and_fit_arbutus2(tree_OUwie, dat1, model = "MV"))
MA_adequacy <- replicate(1000, sim_and_fit_arbutus2(tree_OUwie, dat1, model = "MA"))
MVA_adequacy <- replicate(1000, sim_and_fit_arbutus2(tree_OUwie, dat1, model = "MVA"))

#retrieve pvals
OU_pvals <- OU_adequacy[1,]
MV_pvals <- MV_adequacy[1,]
MA_pvals <- MA_adequacy[1,]
MVA_pvals <- MVA_adequacy[1,]

#arbutus_transform() from custom_functions.R

#transform
OU_df <- arbutus_transform(OU_pvals, 1000) %>% mutate(model = "OU")
MV_df <- arbutus_transform(MV_pvals, 1000) %>% mutate(model = "MV")
MA_df <- arbutus_transform(MA_pvals, 1000) %>% mutate(model = "MA")
MVA_df <- arbutus_transform(MVA_pvals, 1000) %>% mutate(model = "MVA")

#fuse
fuse_df <- full_join(OU_df, MV_df) %>% full_join(MA_df) %>% full_join(MVA_df)

#pivot and plot
data_plot <- fuse_df %>% pivot_longer(cols = c(-model), names_to = "test.stat") %>%
  ggplot(aes(x = value, fill = model)) + geom_boxplot() + facet_wrap(~test.stat)

data_plot

