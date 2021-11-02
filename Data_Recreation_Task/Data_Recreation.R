#Recreating data from MODELING STABILIZING SELECTION: EXPANDING THE ORNSTEIN–UHLENBECK MODEL OF ADAPTIVE EVOLUTION

#First load the required packages
library(geiger)
library(arbutus)
library(ape)
library(tidyverse)
library(OUwie)

#Table 1 data

#First step: Create specific tree topologies. Will be in the phylo format. Comprised of 32, 64, 128, 512 taxa.

#To generate special trees, will use the ape:stree() function. Random tree will be generated using sim.bdtree()

#star tree topology
multistar <- c(stree(32, type = "star"), stree(64, type = "star"), 
               stree(128, type = "star"), stree(512, type = "star"))

##Visualize star topology trees in fan mode because it is more meaningful
starplot <- plot(multistar, type = "fan", layout = 4)

#balanced tree
multibal <- c(stree(32, type = "balanced"), stree(64, type = "balanced"),
              stree(128, type = "balanced"), stree(512, type = "balanced"))
balplot <- plot(multibal, layout = 4)

#pectinate tree will be a left-facing unbalanced tree. Basically the same as right-facing. 
multipec <- c(stree(32, "left"), stree(64, "left"), 
              stree(128, "left"), stree(512, "left"))
pecplot <- plot(multipec, layout = 4)

#random tree will use a birth-death model with b=0.4 and d=0.2
multirand <- c(rphylo(32, 0.4, 0.2), rphylo(64, 0.4, 0.2),
               rphylo(128, 0.4, 0.2), rphylo(512, 0.4, 0.2))
randplot <- plot(multirand, layout = 4)

#Second step: divide selective regime

#For this step, need to make trait dataframes for each phylo object of format 1 species name and 2 selective regime

#star genus and species
star32 <- data.frame(multistar[[1]]$tip.label) %>% rename(Genus_species = multistar..1...tip.label)
star64 <- data.frame(multistar[[2]]$tip.label) %>% rename(Genus_species = multistar..2...tip.label)
star128 <- data.frame(multistar[[3]]$tip.label) %>% rename(Genus_species = multistar..3...tip.label)
star512 <- data.frame(multistar[[4]]$tip.label) %>% rename(Genus_species = multistar..4...tip.label)

#star selective regimes
foo <- c(1:32)
star32 <- star32 %>% 
  mutate(n = foo) %>%
  mutate(Reg = ifelse(n > 16, 2, 1)) %>%
  select(Genus_species, Reg)

foo <- c(1:64)
star64 <- star64 %>% 
  mutate(n = foo) %>%
  mutate(Reg = ifelse(n > 32, 2, 1)) %>%
  select(Genus_species, Reg)

foo <- c(1:128)
star128 <- star128 %>% 
  mutate(n = foo) %>%
  mutate(Reg = ifelse(n > 64, 2, 1)) %>%
  select(Genus_species, Reg)

foo <- c(1:512)
star512 <- star512 %>% 
  mutate(n = foo) %>%
  mutate(Reg = ifelse(n > 256, 2, 1)) %>%
  select(Genus_species, Reg)

#edgelengths
multistar[[1]] <- compute.brlen(multistar[[1]])
multistar[[2]] <- compute.brlen(multistar[[2]])
multistar[[3]] <- compute.brlen(multistar[[3]])
multistar[[4]] <- compute.brlen(multistar[[4]])

multibal[[1]] <- compute.brlen(multibal[[1]])

#Nodes
multistar[[1]] <- makeNodeLabel(multistar[[1]], method = "number", prefix = "")
multistar[[2]] <- makeNodeLabel(multistar[[2]], method = "number", prefix = "")
multistar[[3]] <- makeNodeLabel(multistar[[3]], method = "number", prefix = "")
multistar[[4]] <- makeNodeLabel(multistar[[4]], method = "number", prefix = "")

#Simulate OU models on star topologies
#OUMV
star32_MV <- OUwie.sim(multistar[[1]], star32, alpha = c(1.0,1.0), sigma.sq = c(0.45, 0.90), theta0 = 1.0, theta = c(1.0,2.0))
star64_MV <- OUwie.sim(multistar[[2]], star64, alpha = c(1.0,1.0), sigma.sq = c(0.45, 0.90), theta0 = 1.0, theta = c(1.0,2.0))
star128_MV <- OUwie.sim(multistar[[3]], star128, alpha = c(1.0,1.0), sigma.sq = c(0.45, 0.90), theta0 = 1.0, theta = c(1.0,2.0))
star512_MV <- OUwie.sim(multistar[[4]], star512, alpha = c(1.0,1.0), sigma.sq = c(0.45, 0.90), theta0 = 1.0, theta = c(1.0,2.0))

#OUMA
star32_MA <- OUwie.sim(multistar[[1]], star32, alpha = c(1.0,0.5), sigma.sq = c(0.45, 0.45), theta0 = 1.0, theta = c(1.0,2.0))
star64_MA <- OUwie.sim(multistar[[2]], star64, alpha = c(1.0,0.5), sigma.sq = c(0.45, 0.45), theta0 = 1.0, theta = c(1.0,2.0))
star128_MA <- OUwie.sim(multistar[[3]], star128, alpha = c(1.0,0.5), sigma.sq = c(0.45, 0.45), theta0 = 1.0, theta = c(1.0,2.0))
star512_MA <- OUwie.sim(multistar[[4]], star512, alpha = c(1.0,0.5), sigma.sq = c(0.45, 0.45), theta0 = 1.0, theta = c(1.0,2.0))

#OUMVA
star32_MVA <- OUwie.sim(multistar[[1]], star32, alpha = c(1.0,0.5), sigma.sq = c(0.45, 0.90), theta0 = 1.0, theta = c(1.0,2.0))
star64_MVA <- OUwie.sim(multistar[[2]], star64, alpha = c(1.0,0.5), sigma.sq = c(0.45, 0.90), theta0 = 1.0, theta = c(1.0,2.0))
star128_MVA <- OUwie.sim(multistar[[3]], star128, alpha = c(1.0,0.5), sigma.sq = c(0.45, 0.90), theta0 = 1.0, theta = c(1.0,2.0))
star512_MVA <- OUwie.sim(multistar[[4]], star512, alpha = c(1.0,0.5), sigma.sq = c(0.45, 0.90), theta0 = 1.0, theta = c(1.0,2.0))


#Fit corresponding model
#OUMV OUwie fits
OUwie(phy = multibal[[1]], data = star32_MV, model = "OUMV")
