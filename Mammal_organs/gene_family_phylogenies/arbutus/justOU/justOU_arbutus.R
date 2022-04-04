#Load required libraries
library(ape)
library(geiger)
library(arbutus)
library(tidyverse)
library(flipR)
library(parallel)

### Need to filter out phylogenies that cannot be used

#Load gene trees
final_trees <- readRDS(file = "Mammal_organs/gene_family_phylogenies/genefamilytrees")

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

#Standard Error function
standard_error <- function(x) sd(x) / sqrt(length(x))

#Get SE for each gene per body part
br_SE <- br_dat %>% group_by(hsa) %>% rename(Gene = hsa) %>% summarise(Gene, SE = standard_error(across(hsa.br.M.1:gga.br.F.1)))

#Remove unnecessary data from env
rm(br_dat, amniote_RPKM)

#replace gene IDs with species names
convert <- function(phylo){
  rep <- phylo$tip.label %>%
    str_replace_all("ENSG0..........", "Homo sapiens") %>%
    str_replace_all("...PTRG...........", "Pan troglodytes") %>%
    str_replace_all("...GGOG...........", "Gorilla gorilla") %>%
    str_replace_all("...PPYG...........", "Pongo pygmaeus") %>%
    str_replace_all("...MUSG...........", "Mus musculus") %>%
    str_replace_all("...MMUG...........", "Macaca mulatta") %>%
    str_replace_all("...MODG...........", "Monodelphis domestica") %>%
    str_replace_all("...OANG...........", "Ornithorhynchus anatinus") %>%
    str_replace_all("...GALG...........", "Gallus gallus")
  phylo$tip.label <- rep
  phylo
}

#Now remove entries in data and SE
br_avg_dat <- br_avg_dat %>% filter(Gene %in% names(final_trees))
br_SE <- br_SE %>% filter(Gene %in% names(final_trees))

#Now need to flip tables and properly format
format_expr_data <- function (avgdat) {
  temp <- avgdat %>% pull(Gene)
  avgdat <- avgdat %>% ungroup() %>% select(!Gene)
  dat <- flip(avgdat)
  colnames(dat) <- temp
  res <- dat %>% as.matrix()
  res
}

#Running fitcontinuous

runFC <- function ( dat, SE ){
fitResults <- vector(mode = "list", length = ncol(dat))
for(j in 1:ncol(dat)){
  tdf <- treedata(final_trees[[j]], dat[,j], sort = TRUE)
  phy <- tdf$phy
  data <- tdf$data
  fitOU <- fitContinuous(phy, data, SE[[2]][[j]], model = "OU")
  fitResults[[j]] <- fitOU# %>% arbutus() %>% pvalue_arbutus()
}
fitResults
}

arb_and_pvals <- function(fi){
  res <- fi %>% arbutus() %>% pvalue_arbutus()
  res
}

total_process <- function (dat_list){
  avgdat <- dat_list[[1]]
  part <- dat_list[[2]]
  SE <- dat_list[[3]]
  exp <- format_expr_data(avgdat)
  fit <- runFC(exp, SE)
  fit_name <- paste0("Mammal_organs/gene_family_phylogenies/arbutus/justOU/fit_", part)
  saveRDS(fit, fit_name)
  pvals <- fit %>% map_df(arb_and_pvals) %>% select(!m.sig)
  df_name <- paste0("Mammal_organs/gene_family_phylogenies/arbutus/justOU/pvals_", part)
  saveRDS(pvals, df_name)
  pvals %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("Mammal_organs/gene_family_phylogenies/arbutus/justOU/arbutus_", part, ".png")
  ggsave(pval_name)
}

br_list <- list(br_avg_dat, "br", br_SE)
total_process(br_list)