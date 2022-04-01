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
  bms <- OUwie(tree, df, model = "BMS") %>% arbutus()
  bm <- fitContinuous(tree, df_fix, model = "BM") %>% arbutus()
  ou <- fitContinuous(tree, df_fix, model = "OU") %>% arbutus()
  eb <- fitContinuous(tree, df_fix, model = "EB") %>% arbutus()
  res <- list(bms, bm, ou, eb)
  res
}

test1 <- sim_and_fit_arbutus(tree_OUwie, dat1)

#run the sims
OU_adequacy <- replicate(1000, sim_and_fit_arbutus2(tree_OU, dat2, model = "OU"))
MV_adequacy <- replicate(1000, sim_and_fit_arbutus2(tree_OUwie, dat1, model = "MV"))
MA_adequacy <- replicate(1000, sim_and_fit_arbutus2(tree_OUwie, dat1, model = "MA"))
MVA_adequacy <- replicate(1000, sim_and_fit_arbutus2(tree_OUwie, dat1, model = "MVA"))
BM_adequacy <- replicate(1000, sim_and_fit_arbutus2(tree_OU, dat2, model = "BM"))
EB_adequacy <- replicate(1000, sim_EB(tree_OU, rescaled1, dat2))


#retrieve pvals
OU_pvals <- OU_adequacy[1,]
MV_pvals <- MV_adequacy[1,]
MA_pvals <- MA_adequacy[1,]
MVA_pvals <- MVA_adequacy[1,]
BM_pvals <- BM_adequacy[1,]
EB_pvals <- EB_adequacy[1,]

#arbutus_transform() from custom_functions.R
arbutus_transform <- function ( pval , len , tib) {
  df <- t(data.frame(pval))
  if(missing(tib)) tib = FALSE
  ifelse(tib == TRUE, {
    fin <- as_tibble(df)
  }
  , fin <- data.frame(df))
  row.names(fin) <- c(1:len)
  fin
}

#transform
OU_df <- arbutus_transform(OU_pvals, 1000) %>% mutate(model = "OU")
MV_df <- arbutus_transform(MV_pvals, 1000) %>% mutate(model = "MV")
MA_df <- arbutus_transform(MA_pvals, 1000) %>% mutate(model = "MA")
MVA_df <- arbutus_transform(MVA_pvals, 1000) %>% mutate(model = "MVA")
BM_df <- arbutus_transform(BM_pvals, 1000) %>% mutate(model = "BM")
EB_df <- arbutus_transform(EB_pvals, 1000) %>% mutate(model = "EB")

#fuse
fuse_df <- full_join(OU_df, MV_df) %>% full_join(MA_df) %>% full_join(MVA_df) %>% full_join(BM_df) %>% full_join(EB_df)

#pivot and plot
violin <- fuse_df %>% pivot_longer(cols = c(-model), names_to = "test.stat") %>%
  ggplot(aes(y = value, x = model, fill = (model))) + geom_violin() + geom_boxplot(width = 0.5) + facet_wrap(~test.stat) + theme_bw()

saveRDS(fuse_df,"./Arbutus_Exploration/RDSfiles/Exploration3_data")

ggsave("Arbutus_Exploration/violin_all_models.png", plot = violin)

#reload data
expldata <- readRDS("./Arbutus_Exploration/RDSfiles/Exploration3_data")
#example plot for just cvar
cvarplot <- expldata %>% select(c.var, model) %>% ggplot(aes(x = c.var)) + geom_histogram() + facet_wrap(~model) + theme_bw() 
