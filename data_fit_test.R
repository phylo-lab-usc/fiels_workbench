#Testing fitContinuous with OU models

#First, load required packages

library(tidyverse)
library(ape)
library(geiger)
library(OUwie)
library(broom)

#Now, simulate birth-death tree

tr <- sim.bdtree(n = 128)
plot(tr)

#First test: Simulate character data once, and then fitContinuous 1000 times, with and without fixing alpha bounds

#simulate data use rTraitCont
char_df <- data.frame(rTraitCont(tr, model = "OU"))

#not fixing bounds
OU_fit1 <- replicate(10, fitContinuous(tr, char_df, model = "OU"))  

#fixing bounds
OU_fit2 <- replicate(10, fitContinuous(tr, char_df, model = "OU", bounds = list(alpha = c(0,1.1))))

#What we learned: fitContinuous gets the same result every time, so need to simulate character data

#Second test: Simulate character data 1000 times, then fit continuous to all
sim_and_fit <- function (tree) {
  df <- data.frame(rTraitCont(tree, model = "OU"))
  fit <- fitContinuous(tr, df, model = "OU")
  fit$opt
}

char_list <- replicate(1000, sim_and_fit(tr))
char_df2 <- fix_data_frame(t(data.frame(char_list)))
char_df2$alpha <- unlist(char_df2$alpha)
char_df2 %>% ggplot(aes(x = alpha)) + geom_boxplot() + theme_classic() + annotate("text", x = 2, y = -0.2, label = mean(char_df2$alpha))

#Third test: Simulate using OUwie

dat1 <- data.frame(tr$tip.label) %>% mutate(Reg = 1) %>% rename(Genus_species = tr.tip.label)
tree_OUwie <- tr
tree_OUwie$node.label <- rep(1, 127)


sim_and_fit2 <- function (tree, dat) {
  df <- OUwie.sim(tree, data = dat, alpha = c(1.0, 1.0), sigma.sq = c(0.9, 0.9), theta0 = 1, theta = c(1.0, 1.0))
  df_fix <- df
  row.names(df_fix) <- df_fix$Genus_species
  df_fix <- df_fix %>% select(X)
  fit <- fitContinuous(tree, df_fix, model = "OU")
  fit$opt
}

char_list_v2 <- replicate(1000, sim_and_fit2(tree_OUwie, dat1))
char_df2_v2 <- fix_data_frame(t(data.frame(char_list_v2)))
char_df2_v2$alpha <- unlist(char_df2_v2$alpha)

char_df2_v2 %>% ggplot(aes(x = alpha)) + geom_boxplot() + theme_classic() + annotate("text", x = 2, y = -0.2, label = mean(char_df2$alpha))
