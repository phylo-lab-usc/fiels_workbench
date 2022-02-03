#Testing fitContinuous with OU models

#First, load required packages

library(tidyverse)
library(ape)
library(geiger)
library(OUwie)
library(broom)
library(arbutus)

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


#Run adequacy tests on fits, this should be the standard to compare to other OU models
#previously only got the opt data from sims, need the whole fit objects

#Fit 1
sim_and_fit_arbutus <- function (tree) {
  df <- data.frame(rTraitCont(tree, model = "OU"))
  fit <- fitContinuous(tr, df, model = "OU")
  a <- arbutus(fit)
  a
}



#Fit 2
sim_and_fit_arbutus2 <- function (tree, dat) {
  df <- OUwie.sim(tree, data = dat, alpha = c(1.0, 1.0), sigma.sq = c(0.9, 0.9), theta0 = 1, theta = c(1.0, 1.0))
  df_fix <- df
  row.names(df_fix) <- df_fix$Genus_species
  df_fix <- df_fix %>% select(X)
  fit <- fitContinuous(tree, df_fix, model = "OU")
  a <- arbutus(fit)
  a
}

#run the sims
first_sim_adequacy <- replicate(1000, sim_and_fit_arbutus(tr))
second_sim_adequacy <- replicate(1000, sim_and_fit_arbutus2(tr, dat1))

#data transform function
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

#retrieve p values
first_sim_pvals <- first_sim_adequacy[1,]
first_sim_dfs <- arbutus_transform(first_sim_pvals, 1000)
second_sim_dfs <- arbutus_transform(second_sim_adequacy[1,], 1000)

#edit data to plot

#first label each df
first_sim_dfs_fuse <- first_sim_dfs %>% mutate(lab = "first")
second_sim_dfs_fuse <- second_sim_dfs %>% mutate(lab = "second")

#Now stack on top
sims_df <- full_join(first_sim_dfs_fuse, second_sim_dfs_fuse)

#Now pivot
sims_df_piv <- sims_df %>% 
  pivot_longer(cols = c(-lab), names_to = "test.stat")

#now plot
sims_df_piv %>% group_by(test.stat) %>% ggplot(aes(x = value, fill = lab)) + geom_boxplot() + facet_wrap(~test.stat)

#shows that OUwie sim and rTraitCont work equally the same. Can now use this method with other models
