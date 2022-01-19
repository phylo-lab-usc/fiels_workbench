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
cb_dat <- amniote_RPKM %>% select(hsa,contains(".cb."))
ht_dat <- amniote_RPKM %>% select(hsa,contains(".ht."))
kd_dat <- amniote_RPKM %>% select(hsa,contains(".kd."))
lv_dat <- amniote_RPKM %>% select(hsa,contains(".lv."))
ts_dat <- amniote_RPKM %>% select(hsa,contains(".ts."))

#Take averages of each species for each body [art]
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

cb_avg_dat <- cb_dat %>% group_by(hsa) %>% rename(Gene = hsa) %>% transmute("Homo sapiens" = mean(hsa.cb.M.1, hsa.cb.F.1),
                                                                            "Pan troglodytes" = mean(ptr.cb.M.1, ptr.cb.F.1),
                                                                            "Pan paniscus" = mean(ppa.cb.M.1, ppa.cb.F.1),
                                                                            "Gorilla gorilla" = mean(ggo.cb.M.1, ggo.cb.F.1),
                                                                            "Pongo pygmaeus" = mean(ppy.cb.F.1),
                                                                            "Mus musculus" = mean(mmu.cb.M.1, mmu.cb.M.2, mmu.cb.F.1),
                                                                            "Macaca mulatta" = mean(mml.cb.F.1, mml.cb.M.1),
                                                                            "Monodelphis domestica" = mean(mdo.cb.M.1, mdo.cb.F.1),
                                                                            "Ornithorhynchus anatinus" = mean(oan.cb.M.1, oan.cb.F.1),
                                                                            "Gallus gallus" = mean(gga.cb.M.1, gga.cb.F.1))

ht_avg_dat <- ht_dat %>% group_by(hsa) %>% rename(Gene = hsa) %>% transmute("Homo sapiens" = mean(hsa.ht.M.1, hsa.ht.M.2, hsa.ht.F.1),
                                                                            "Pan troglodytes" = mean(ptr.ht.M.1, ptr.ht.F.1),
                                                                            "Pan paniscus" = mean(ppa.ht.M.1, ppa.ht.F.1),
                                                                            "Gorilla gorilla" = mean(ggo.ht.M.1, ggo.ht.F.1),
                                                                            "Pongo pygmaeus" = mean(ppy.ht.M.1, ppy.ht.F.1),
                                                                            "Mus musculus" = mean(mmu.ht.M.1, mmu.ht.M.2, mmu.ht.F.1),
                                                                            "Macaca mulatta" = mean(mml.ht.F.1, mml.ht.M.1),
                                                                            "Monodelphis domestica" = mean(mdo.ht.M.1, mdo.ht.F.1),
                                                                            "Ornithorhynchus anatinus" = mean(oan.ht.M.1, oan.ht.F.1),
                                                                            "Gallus gallus" = mean(gga.ht.M.1, gga.ht.F.1))

kd_avg_dat <- kd_dat %>% group_by(hsa) %>% rename(Gene = hsa) %>% transmute("Homo sapiens" = mean(hsa.kd.M.1, hsa.kd.M.2, hsa.kd.F.1),
                                                                            "Pan troglodytes" = mean(ptr.kd.M.1, ptr.kd.F.1),
                                                                            "Pan paniscus" = mean(ppa.kd.M.1, ppa.kd.M.1),
                                                                            "Gorilla gorilla" = mean(ggo.kd.M.1, ggo.kd.F.1),
                                                                            "Pongo pygmaeus" = mean(ppy.kd.M.1, ppy.kd.F.1),
                                                                            "Mus musculus" = mean(mmu.kd.M.1, mmu.kd.M.2, mmu.kd.F.1),
                                                                            "Macaca mulatta" = mean(mml.kd.F.1, mml.kd.M.1),
                                                                            "Monodelphis domestica" = mean(mdo.kd.M.1, mdo.kd.F.1),
                                                                            "Ornithorhynchus anatinus" = mean(oan.kd.M.1, oan.kd.F.1),
                                                                            "Gallus gallus" = mean(gga.kd.M.1, gga.kd.F.1))

lv_avg_dat <- lv_dat %>% group_by(hsa) %>% rename(Gene = hsa) %>% transmute("Homo sapiens" = mean(hsa.lv.M.1, hsa.lv.M.2),
                                                                            "Pan troglodytes" = mean(ptr.lv.M.1, ptr.lv.F.1),
                                                                            "Pan paniscus" = mean(ppa.lv.F.1, ppa.lv.F.1),
                                                                            "Gorilla gorilla" = mean(ggo.lv.M.1, ggo.lv.F.1),
                                                                            "Pongo pygmaeus" = mean(ppy.lv.M.1, ppy.lv.F.1),
                                                                            "Mus musculus" = mean(mmu.lv.M.1, mmu.lv.M.2, mmu.lv.F.1),
                                                                            "Macaca mulatta" = mean(mml.lv.F.1, mml.lv.M.1),
                                                                            "Monodelphis domestica" = mean(mdo.lv.M.1, mdo.lv.F.1),
                                                                            "Ornithorhynchus anatinus" = mean(oan.lv.M.1, oan.lv.F.1),
                                                                            "Gallus gallus" = mean(gga.lv.M.1, gga.lv.F.1))

ts_avg_dat <- ts_dat %>% group_by(hsa) %>% rename(Gene = hsa) %>% transmute("Homo sapiens" = mean(hsa.ts.M.1, hsa.ts.M.2),
                                                                            "Pan troglodytes" = mean(ptr.ts.M.1),
                                                                            "Pan paniscus" = mean(ppa.ts.M.1),
                                                                            "Gorilla gorilla" = mean(ggo.ts.M.1),
                                                                            "Mus musculus" = mean(mmu.ts.M.1, mmu.ts.M.2),
                                                                            "Macaca mulatta" = mean(mml.ts.M.1, mml.ts.M.2),
                                                                            "Monodelphis domestica" = mean(mdo.ts.M.1, mdo.ts.M.2),
                                                                            "Ornithorhynchus anatinus" = mean(oan.ts.M.1, oan.ts.M.2, oan.ts.M.3),
                                                                            "Gallus gallus" = mean(gga.ts.M.1, gga.ts.M.2))

#Remove unnecessary data from env
rm(br_dat, cb_dat, ht_dat, kd_dat, lv_dat, ts_dat, amniote_RPKM)

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
  fitBM <- fitContinuous(phy, data[,1], model = "BM")
  fitOU <- fitContinuous(phy, data[,1], model = "OU")
  fitEB <- fitContinuous(phy, data[,1], model = "EB")
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
  #fit_name <- paste0("Mammal_organs/species_phylogeny/arbutus/fit_", part)
  #saveRDS(fit, file = fit_name)
  df <- model_count(fit)
  aic_name <- paste0("Mammal_organs/species_phylogeny/AIC/AIC_", part, ".png")
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  ggsave(aic_name)
  result <- run_arb(fit)
  rds_name <- paste0("Mammal_organs/species_phylogeny/arbutus/pvals_", part)
  saveRDS(result, file = rds_name)
  result %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("Mammal_organs/species_phylogeny/arbutus/arbutus_", part, ".png")
  ggsave(pval_name)
}

total_process(br_avg_dat, "br")
total_process(br_avg_dat, "cb")
total_process(br_avg_dat, "ht")
total_process(br_avg_dat, "kd")
total_process(br_avg_dat, "lv")
total_process(br_avg_dat, "ts")