#Test script

library(tidyverse)
library(geiger)
library(arbutus)
library(ape)
library(OUwie)

#data(finch)
#geo <- finch$phy
#dat <- finch$data[,"wingL"]

#first_sim <- sim.char(geo, 0.02, 1000)

#unit_tree <- make_unit_tree(geo, data = dat)

#first_sim_unit <- simulate_char_unit(unit_tree)

#test.stat <- calculate_pic_stat(unit_tree)


#Testing how to scale Early Burst trees. 
test.tree <- sim.bdtree(n = 128)
data <- data.frame(test.tree$tip.label) %>% mutate(Reg = 1) %>% rename(Genus_species = test.tree.tip.label)
test.tree$node.label <- rep(1, 127)

#testing parameter a
rescaled1 <- rescale(test.tree, model = "EB", a = 1)
rescaledhalf <- rescale(test.tree, model = "EB", a = 0.5)
rescaled2 <- rescale(test.tree, model = "EB", a = 2)

sim_and_fit_EB <- function (tree, dat) {
  df <- OUwie.sim(tree, dat, alpha = c(1e-10, 1e-10), sigma.sq = c(0.45, 0.45), theta0 = 1.0, theta = c(0, 0))
  df_fix <- df
  row.names(df_fix) <- df_fix$Genus_species
  df_fix <- df_fix %>% select(X)
  fit <- fitContinuous(tree, df_fix, model = "OU")
  a <- arbutus(fit)
  a
}

#run simulations
BMhalf <- replicate(1000, sim_and_fit_EB(rescaledhalf, data))
BM1 <- replicate(1000, sim_and_fit_EB(rescaled1, data))
BM2 <- replicate(1000, sim_and_fit_EB(rescaled2, data))

#retrieve pvals
BMhalf_pvals <- BMhalf[1,]
BM1_pvals <- BM1[1,]
BM2_pvals <- BM2[1,]


#arbutus_transform() from custom_functions.R

#transform
BMhalf_df <- arbutus_transform(BMhalf_pvals, 1000) %>% mutate(alpha = "half")
BM1_df <- arbutus_transform(BM1_pvals, 1000) %>% mutate(alpha = "one")
BM2_df <- arbutus_transform(BM2_pvals, 1000) %>% mutate(alpha = "two")

#fuse
fuse_df <- full_join(BMhalf_df, BM1_df) %>% full_join(BM2_df) 

#pivot and plot
fuse_df %>% pivot_longer(cols = c(-alpha), names_to = "test.stat") %>%
  ggplot(aes(y = value, x = alpha, fill = (alpha))) + geom_violin() + geom_boxplot(width = 0.5) + facet_wrap(~test.stat) + theme_bw()
