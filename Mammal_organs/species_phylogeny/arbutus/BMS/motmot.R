#Load libraries
library(motmot)
library(tidyverse)
library(geiger)
library(flipR)

#Load data
species_phylo <- read.tree(file = "Mammal_organs/species_phylogeny/species_names.nwk")

#Rename tips
species_phylo$tip.label <- gsub("_", " ", species_phylo$tip.label)

#Load gene expression data
amniote_RPKM <- read.delim("Mammal_organs/Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.txt")

#Split into body parts
br_dat <- amniote_RPKM %>% select(hsa,contains(".br.")) 

#Take averages of each species for each body part
br_avg_dat <- br_dat %>% group_by(hsa) %>% rename(Gene = hsa) %>% transmute("Homo sapiens" = mean(hsa.br.M.1, hsa.br.M.2, hsa.br.M.3, hsa.br.M.4, hsa.br.M.5, hsa.br.F.1),
                                                                            "Pan troglodytes" = mean(ptr.br.M.1, ptr.br.M.2, ptr.br.M.3, ptr.br.M.4, ptr.br.M.5, ptr.br.F.1),
                                                                            "Pan paniscus" = mean(ppa.br.M.1, ppa.br.F.2, ppa.br.F.1),
                                                                            "Gorilla gorilla" = mean(ggo.br.M.1, ggo.br.F.1),
                                                                            "Pongo pygmaeus" = mean(ppy.br.M.1, ppy.br.F.1),
                                                                            "Mus musculus" = mean(mmu.br.M.1, mmu.br.M.2, mmu.br.F.1),
                                                                            "Macaca mulatta" = mean(mml.br.F.1, mml.br.M.1, mml.br.M.2),
                                                                            "Monodelphis domestica" = mean(mdo.br.M.1, mdo.br.F.1),
                                                                            "Ornithorhynchus anatinus" = mean(oan.br.M.1, oan.br.F.1),
                                                                            "Gallus gallus" = mean(gga.br.M.1, gga.br.F.1))
#Now need to flip tables and properly format
format_expr_data <- function (avgdat) {
  temp <- avgdat %>% pull(Gene)
  avgdat <- avgdat %>% ungroup() %>% select(!Gene)
  dat <- flip(avgdat)
  colnames(dat) <- temp
  dat
}

data <- br_avg_dat %>% format_expr_data() %>% as.matrix()
#Now split the brain data into pieces
motmot_df <- data[,1]

td <- treedata(species_phylo, motmot_df)

#Find rate shifts
tm2 <- transformPhylo.ML(td$data, td$phy, model = "tm2", nSplits = 3)
sum_tm2 <- summary(tm2)
plot(sum_tm2)
cutOff <- calcCutOff(td$phy, n = 1000, mc.cores = 4, model = "tm2")
res1 <- summary(tm2, cutoff = cutOff)
saveRDS(res1, "Mammal_organs/species_phylogeny/arbutus/BMS/motmot.rds")