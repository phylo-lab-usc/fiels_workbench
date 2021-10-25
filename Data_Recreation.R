#Recreating data from MODELING STABILIZING SELECTION: EXPANDING THE ORNSTEIN–UHLENBECK MODEL OF ADAPTIVE EVOLUTION

#First load the required packages
library(geiger)
library(arbutus)
library(ape)

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

#star trees were split between regime 1 and 2
