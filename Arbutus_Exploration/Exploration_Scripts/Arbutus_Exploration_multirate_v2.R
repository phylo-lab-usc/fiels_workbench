#Multirate arbutus test 2: True model is BM, testing if AIC will detect BM vs BMS

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
  df <- OUwie.sim(tree, dat2, alpha = c(1e-10, 1e-10), sigma.sq = c(0.45, 0.45), theta0 = 1.0, theta = c(0, 0))
  df_fix <- df
  row.names(df_fix) <- df_fix$Genus_species
  df_fix <- df_fix %>% select(X)
  bms <- OUwie(tree, df, model = "BMS") 
  bm <- fitContinuous(tree, df_fix, model = "BM") 
  ou <- fitContinuous(tree, df_fix, model = "OU") 
  eb <- fitContinuous(tree, df_fix, model = "EB") 
  aic <- c(bm$opt[["aic"]], ou$opt[["aic"]], eb$opt[["aic"]], bms$AIC)
  res <- ifelse(min(aic) == aic[1], list(c(bm, model = "BM")), 
                ifelse(min(aic) == aic[2], list(c(ou, model = "OU")), 
                       ifelse(min(aic) == aic[3], list(c(eb, model = "EB")),
                              list(c(bms, model = "BMS")))))
  res
}

#run the sims
AICs <- replicate(1000, sim_and_fit_arbutus(tree_OUwie, dat1))
saveRDS(AICs, "Arbutus_Exploration/RDSfiles/AIC_multirate.rds")

model_count <- function (fit) {
  ou = 0
  bm = 0
  eb = 0
  bms = 0
  for(f in fit){
    vec <- f
    ifelse(vec$model == "OU", ou <- ou + 1, ifelse(vec$model == "BM", bm <- bm + 1, ifelse(vec$model == "EB", eb <- eb + 1, bms <- bms + 1)))
  }
  df <- data.frame(OU = ou, BM = bm, EB = eb, BMS = bms)
  b <- df %>% pivot_longer(c(OU, BM, EB, BMS), names_to = "model")
  b
}

model_count(AICs) %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
ggsave("Arbutus_Exploration/Figures/Multirate_AIC.png")
