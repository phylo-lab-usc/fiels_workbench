#Arbutus exploration via OU simulation

#Task: Simulate OU models (MV, MA, MVA) and then fit simple OU model to them and check adequacy.

#load required packages
library(ape)
library(geiger)
library(arbutus)
library(OUwie)
library(tidyverse)

#First, make tree topology. For first analysis will use same tree for all. 

tr <- sim.bdtree(n = 128) %>% 
  compute.brlen()

#regime labels needed for next step
regime <- c(rep(1, 64), rep(2, 64))
tr$node.label <- regime[1:(length(regime) -1)]

#Now create data frame

df <- data.frame(tr$tip.label) %>% 
  rename(Genus_species = tr.tip.label) %>%
  mutate(Reg = regime)

#Now for simulation of each OU model
df_MV <- OUwie.sim(tr, df, alpha = c(1.0, 1.0), sigma.sq = c(0.45, 0.9), theta0 = 1.0, theta = c(1.0, 2.0))
df_MA <- OUwie.sim(tr, df, alpha = c(1.0, 0.5), sigma.sq = c(0.45, 0.45), theta0 = 1.0, theta = c(1.0, 2.0))
df_MVA <- OUwie.sim(tr, df, alpha = c(1.0, 0.5), sigma.sq = c(0.45, 0.9), theta0 = 1.0, theta = c(1.0, 2.0))

#Now to pivot data for use in fitting, need to fit data according to dat format, pull format from below
geo <- get(data("geospiza"))
geo$dat

#Fixing dfs

df_MV_fix <- df_MV
row.names(df_MV_fix) <- df_MV$Genus_species
df_MV_fix <- df_MV_fix %>% select(X)

df_MA_fix <- df_MA
row.names(df_MA_fix) <- df_MA$Genus_species
df_MA_fix <- df_MA_fix %>% select(X)

df_MVA_fix <- df_MVA
row.names(df_MVA_fix) <- df_MVA$Genus_species
df_MVA_fix <- df_MVA_fix %>% select(X)

#Now fit a simple OU model to each OU type

fit_MV <- fitContinuous(tr, df_MV_fix, model = "OU")
fit_MA <- fitContinuous(tr, df_MA_fix, model = "OU")
fit_MVA <- fitContinuous(tr, df_MVA_fix, model = "OU")

#Now for adequacy

#First make unit trees
unit_MV <- make_unit_tree(fit_MV)
unit_MA <- make_unit_tree(fit_MA)
unit_MVA <- make_unit_tree(fit_MVA)

#Now run arbutus
arby <- function(unit.tree) {
  calc1 <- calculate_pic_stat(unit.tree)
  sim.data <- simulate_char_unit(unit.tree)
  calc2 <- calculate_pic_stat(sim.data)
  compare_pic_stat(calc1, calc2)
}

#for comparison, will simulate a regular OU model and then calculate adequacy
df_OU <- data.frame(rTraitCont(tr, model = "OU"))
fit_OU <- fitContinuous(tr, df_OU, model = "OU")
unit_OU <- make_unit_tree(fit_OU)
OU_adequacy <- arby(unit_OU)

#Adequacy simulations
adequacy <- function(){
  MV_adequacy <- arby(unit_MV)
  MA_adequacy <- arby(unit_MA)
  MVA_adequacy <- arby(unit_MVA)
  OU_adequacy <- arby(unit_OU)
  list(MV_adequacy$p.values, MA_adequacy$p.values, MVA_adequacy$p.values, OU_adequacy$p.values)
}

overall <- replicate(100, adequacy())
