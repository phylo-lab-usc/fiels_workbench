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
amniote_RPKM <- read.delim("Mammal_organs/Supplementary_Data1/NormalizedRPKM_ConstitutiveAlignedExons_Amniote1to1Orthologues.txt") %>%
  rename(Gene = hsa)

#First split data into body part, and then species
df_format <- function( body_part ) {
  hsa <- amniote_RPKM %>% select(Gene,contains(body_part)) %>% select(Gene,contains("hsa")) %>%
    pivot_longer(!Gene, names_to = "replicate", values_to = "Homo sapiens") %>% select(!replicate)
  ptr <- amniote_RPKM %>% select(Gene,contains(body_part)) %>% select(Gene,contains("ptr")) %>%
    pivot_longer(!Gene, names_to = "replicate", values_to = "Pan troglodytes") %>% select(!replicate)
  ppa <- amniote_RPKM %>% select(Gene,contains(body_part)) %>% select(Gene,contains("ppa")) %>%
    pivot_longer(!Gene, names_to = "replicate", values_to = "Pan paniscus") %>% select(!replicate)
  ggo <- amniote_RPKM %>% select(Gene,contains(body_part)) %>% select(Gene,contains("ggo")) %>%
    pivot_longer(!Gene, names_to = "replicate", values_to = "Gorilla gorilla") %>% select(!replicate)
  ppy <- amniote_RPKM %>% select(Gene,contains(body_part)) %>% select(Gene,contains("ppy")) %>%
    pivot_longer(!Gene, names_to = "replicate", values_to = "Pongo pygmaeus") %>% select(!replicate)
  mml <- amniote_RPKM %>% select(Gene,contains(body_part)) %>% select(Gene,contains("mml")) %>%
    pivot_longer(!Gene, names_to = "replicate", values_to = "Macaca mulatta") %>% select(!replicate)
  mmu <- amniote_RPKM %>% select(Gene,contains(body_part)) %>% select(Gene,contains("mmu")) %>%
    pivot_longer(!Gene, names_to = "replicate", values_to = "Mus musculus") %>% select(!replicate)
  mdo <- amniote_RPKM %>% select(Gene,contains(body_part)) %>% select(Gene,contains("mdo")) %>%
    pivot_longer(!Gene, names_to = "replicate", values_to = "Monodelphis domestica") %>% select(!replicate)
  oan <- amniote_RPKM %>% select(Gene,contains(body_part)) %>% select(Gene,contains("oan")) %>%
    pivot_longer(!Gene, names_to = "replicate", values_to = "Ornithorhynchus anatinus") %>% select(!replicate)
  gga <- amniote_RPKM %>% select(Gene,contains(body_part)) %>% select(Gene,contains("gga")) %>%
    pivot_longer(!Gene, names_to = "replicate", values_to = "Gallus gallus") %>% select(!replicate)
  
  final <- inner_join(hsa, ptr) %>% inner_join(ppa) %>% inner_join(ggo) %>% inner_join(ppy) %>%
    inner_join(mml) %>% inner_join(mmu) %>% inner_join(mdo) %>% inner_join(oan) %>% inner_join(gga)
  final
}

#Split into body parts
br_dat <- df_format( "br" )
cb_dat <- df_format( "cb" )
ht_dat <- df_format( "ht" )
kd_dat <- df_format( "kd" )
lv_dat <- df_format( "lv" )
ts_dat <- df_format( "ts" )

#Remove unnecessary data from env
rm(amniote_RPKM)

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
  df <- model_count(fit)
  aic_name <- paste0("Mammal_organs/species_phylogeny/AIC/non_avg_dat/AIC_", part, ".png")
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  ggsave(aic_name)
  result <- run_arb(fit)
  rds_name <- paste0("Mammal_organs/species_phylogeny/arbutus/non_avg_dat/pvals_", part)
  saveRDS(result, file = rds_name)
  result %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("Mammal_organs/species_phylogeny/arbutus/non_avg_dat/arbutus_", part, ".png")
  ggsave(pval_name)
}

total_process(br_dat, "br")
total_process(cb_dat, "cb")
total_process(ht_dat, "ht")
total_process(kd_dat, "kd")
total_process(lv_dat, "lv")
total_process(ts_dat, "ts")