#Load required libraries
library(ape)
library(geiger)
library(arbutus)
library(tidyverse)
library(flipR)

#Load newick format phylogeny and convert to something recognizable by ape/geiger
species_phylo <- read.tree(file = "Mammal_organs/species_phylogeny/species_names.nwk")

#Rename tips
species_phylo$tip.label <- gsub("_", " ", species_phylo$tip.label)

#Load gene expression data
amniote_RPKM <- read.delim("Mammal_organs/Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.txt")

#Split into body parts
br_dat <- amniote_RPKM %>% select(hsa,contains(".br.")) 

#Take averages of each species for each body part
br_avg_dat <- br_dat %>% rename(Gene = hsa) %>% group_by(Gene) %>% transmute("Homo sapiens" = mean(c(hsa.br.M.1, hsa.br.M.2, hsa.br.M.3, hsa.br.M.4, hsa.br.M.5, hsa.br.F.1)),
                                                  "Pan troglodytes" = mean(c(ptr.br.M.1, ptr.br.M.2, ptr.br.M.3, ptr.br.M.4, ptr.br.M.5, ptr.br.F.1)),
                                                  "Pan paniscus" = mean(c(ppa.br.M.1, ppa.br.F.2, ppa.br.F.1)),
                                                  "Gorilla gorilla" = mean(c(ggo.br.M.1, ggo.br.F.1)),
                                                  "Pongo pygmaeus" = mean(c(ppy.br.M.1, ppy.br.F.1)),
                                                  "Mus musculus" = mean(c(mmu.br.M.1, mmu.br.M.2, mmu.br.F.1)),
                                                  "Macaca mulatta" = mean(c(mml.br.F.1, mml.br.M.1, mml.br.M.2)),
                                                  "Monodelphis domestica" = mean(c(mdo.br.M.1, mdo.br.F.1)),
                                                  "Ornithorhynchus anatinus" = mean(c(oan.br.M.1, oan.br.F.1)),
                                                  "Gallus gallus" = mean(c(gga.br.M.1, gga.br.F.1))) %>% 
  ungroup() %>% slice_sample(n = 100) %>% arrange(Gene)

std.error <- function ( x ) {
  sd(x) / sqrt(length(x))
}

br_SE <- br_dat %>% rename(Gene = hsa) %>% filter(Gene %in% br_avg_dat$Gene) %>% group_by(Gene) %>%
  transmute("Homo sapiens" = std.error(c(hsa.br.M.1, hsa.br.M.2, hsa.br.M.3, hsa.br.M.4, hsa.br.M.5, hsa.br.F.1)),
            "Pan troglodytes" = std.error(c(ptr.br.M.1, ptr.br.M.2, ptr.br.M.3, ptr.br.M.4, ptr.br.M.5, ptr.br.F.1)),
            "Pan paniscus" = std.error(c(ppa.br.M.1, ppa.br.F.2, ppa.br.F.1)),
            "Gorilla gorilla" = std.error(c(ggo.br.M.1, ggo.br.F.1)),
            "Pongo pygmaeus" = std.error(c(ppy.br.M.1, ppy.br.F.1)),
            "Mus musculus" = std.error(c(mmu.br.M.1, mmu.br.M.2, mmu.br.F.1)),
            "Macaca mulatta" = std.error(c(mml.br.F.1, mml.br.M.1, mml.br.M.2)),
            "Monodelphis domestica" = std.error(c(mdo.br.M.1, mdo.br.F.1)),
            "Ornithorhynchus anatinus" = std.error(c(oan.br.M.1, oan.br.F.1)),
            "Gallus gallus" = std.error(c(gga.br.M.1, gga.br.F.1))) %>% ungroup() %>% arrange(Gene)



#Remove unnecessary data from env
rm(br_dat, amniote_RPKM)

#Now need to flip tables and properly format
format_expr_data <- function (avgdat) {
  temp <- avgdat %>% pull(Gene)
  avgdat <- avgdat %>% ungroup() %>% select(!Gene)
  dat <- flip(avgdat)
  colnames(dat) <- temp
  dat
}

#Running fitcontinuous

runFC <- function (dat){
fitResults <- vector(mode = "list", length = ncol(dat))
tdf <- treedata(species_phylo, dat, sort = TRUE)
phy <- tdf$phy
data <- tdf$data
for(j in 1:ncol(dat)){
  fitBM <- fitContinuous(phy, data[,j], model = "BM")
  fitOU <- fitContinuous(phy, data[,j], model = "OU")
  fitEB <- fitContinuous(phy, data[,j], model = "EB")
  aic <- c(fitBM$opt[["aic"]], fitOU$opt[["aic"]], fitEB$opt[["aic"]])
  fit <- ifelse(min(aic) == aic[1], list(c(fitBM, model = "BM")), 
                ifelse(min(aic) == aic[2], list(c(fitOU, model = "OU")), 
                       list(c(fitEB, model = "EB"))))
  fitResults[j] <- fit
}
fitResults
}

model_count <- function (fit) {
  ou = 0
  bm = 0
  eb = 0
  for(f in fit){
    vec <- f
    ifelse(vec$model == "OU", ou <- ou + 1, ifelse(vec$model == "BM", bm <- bm + 1, eb <- eb + 1))
  }
  df <- data.frame(OU = ou, BM = bm, EB = eb)
  b <- df %>% pivot_longer(c(OU, BM, EB), names_to = "model")
  b
}

#running arbutus
run_arb <- function (fits){
  arby <- vector("list", length = length(fits))
  count = 1
  for(f in fits){
    class(f) <- "gfit"
    arby[[count]] <- arbutus(f)
    count = count + 1
  }
  arby_df <- map_df(arby, pvalue_arbutus)
  arby_df
}

total_process <- function (avgdat, part){
  exp <- format_expr_data(avgdat)
  fit <- runFC(exp)
  fit_name <- paste0("Mammal_organs/species_phylogeny/arbutus/fit_", part)
  saveRDS(fit, file = fit_name)
  df <- model_count(fit)
  aic_name <- paste0("Mammal_organs/species_phylogeny/AIC/AIC_", part, ".png")
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  #ggsave(aic_name)
  result <- run_arb(fit)
  rds_name <- paste0("Mammal_organs/species_phylogeny/arbutus/pvals_", part)
  saveRDS(result, file = rds_name)
  result %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("Mammal_organs/species_phylogeny/arbutus/arbutus_", part, ".png")
  #ggsave(pval_name)
}

total_process(br_avg_dat, "br")
