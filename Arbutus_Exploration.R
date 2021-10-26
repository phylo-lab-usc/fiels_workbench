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