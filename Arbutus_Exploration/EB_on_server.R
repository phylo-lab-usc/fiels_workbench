#Test EB fit compared to BM and OU

library(OUwie)
library(tidyverse)
library(here)

#Make tree and data
test.tree <- sim.bdtree(n = 128)
data <- data.frame(test.tree$tip.label) %>% mutate(Reg = 1) %>% rename(Genus_species = test.tree.tip.label)
test.tree$node.label <- rep(1, 127)

#testing parameter a
rescaled1 <- rescale(test.tree, model = "EB", a = 1)
rescaledhalf <- rescale(test.tree, model = "EB", a = 0.5)
rescaled2 <- rescale(test.tree, model = "EB", a = 2)

#Sim and fit to OU
sim_OU <- function (tree, rescaled, dat) {
df <- OUwie.sim(rescaled, dat, alpha = c(1e-10, 1e-10), sigma.sq = c(0.45, 0.45), theta0 = 1.0, theta = c(0, 0))
df_fix <- df
row.names(df_fix) <- df_fix$Genus_species
df_fix <- df_fix %>% select(X)
fit <- fitContinuous(tree, df_fix, model = "OU")
fit$opt
}

#Sim and fit to BM
sim_BM <- function (tree, rescaled, dat) {
  df <- OUwie.sim(rescaled, dat, alpha = c(1e-10, 1e-10), sigma.sq = c(0.45, 0.45), theta0 = 1.0, theta = c(0, 0))
  df_fix <- df
  row.names(df_fix) <- df_fix$Genus_species
  df_fix <- df_fix %>% select(X)
  fit <- fitContinuous(tree, df_fix, model = "BM")
  fit$opt
}

#Sim and fit to EB
sim_EB <- function (tree, rescaled, dat) {
  df <- OUwie.sim(rescaled, dat, alpha = c(1e-10, 1e-10), sigma.sq = c(0.45, 0.45), theta0 = 1.0, theta = c(0, 0))
  df_fix <- df
  row.names(df_fix) <- df_fix$Genus_species
  df_fix <- df_fix %>% select(X)
  fit <- fitContinuous(tree, df_fix, model = "EB")
  fit$opt
}

#Run the simulations fit to OU
OU_ctrl <- replicate(500, sim_OU(test.tree, test.tree, data))
OU_half <- replicate(500, sim_OU(test.tree, rescaledhalf, data))
OU_1 <- replicate(500, sim_OU(test.tree, rescaled1, data))
OU_2 <- replicate(500, sim_OU(test.tree, rescaled2, data))

#Run the simulations fit to BM
BM_ctrl <- replicate(500, sim_BM(test.tree, test.tree, data))
BM_half <- replicate(500, sim_BM(test.tree, rescaledhalf, data))
BM_1 <- replicate(500, sim_BM(test.tree, rescaled1, data))
BM_2 <- replicate(500, sim_BM(test.tree, rescaled2, data))

#Run the simulations fit to EB
EB_half <- replicate(500, sim_EB(test.tree, rescaledhalf, data))
EB_1 <- replicate(500, sim_EB(test.tree, rescaled1, data))
EB_2 <- replicate(500, sim_EB(test.tree, rescaled2, data))

#Function to fix data
fit_transform <- function ( fit , len , tib) {
  df <- t(data.frame(fit[1,]))
  if(missing(tib)) tib = FALSE
  ifelse(tib == TRUE, {
    fin <- as_tibble(df)
  }
  , fin <- data.frame(df))
  row.names(fin) <- c(1:len)
  fin
}


#Get sigma.sq vals
OU_ctrl_sig <- fit_transform(OU_ctrl, 500) %>% mutate(alpha = "ctrl")
OU_half_sig <- fit_transform(OU_half, 500) %>% mutate(alpha = "half")
OU_1_sig <- fit_transform(OU_1, 500) %>% mutate(alpha = "1")
OU_2_sig <- fit_transform(OU_2, 500) %>% mutate(alpha = "2")

BM_ctrl_sig <- fit_transform(BM_ctrl, 500) %>% mutate(alpha = "ctrl")
BM_half_sig <- fit_transform(BM_half, 500) %>% mutate(alpha = "half")
BM_1_sig <- fit_transform(BM_1, 500) %>% mutate(alpha = "1")
BM_2_sig <- fit_transform(BM_2, 500) %>% mutate(alpha = "2")

EB_half_sig <- fit_transform(EB_half, 500) %>% mutate(alpha = "half")
EB_1_sig <- fit_transform(EB_1, 500) %>% mutate(alpha = "1")
EB_2_sig <- fit_transform(EB_2, 500) %>% mutate(alpha = "2")

#fuse
OU_df <- full_join(OU_ctrl_sig, OU_half_sig) %>% full_join(OU_1_sig) %>% full_join(OU_2_sig)
BM_df <- full_join(BM_ctrl_sig, BM_half_sig) %>% full_join(BM_1_sig) %>% full_join(BM_2_sig)
EB_df <- EB_half_sig %>% full_join(EB_1_sig) %>% full_join(EB_2_sig)

#save data for later
saveRDS(OU_df, paste0(here(), "EB_fit_to_OU"))
saveRDS(BM_df, paste0(here(), "EB_fit_to_BM"))
saveRDS(EB_df, paste0(here(), "EB_fit_to_EB"))
