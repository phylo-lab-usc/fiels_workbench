# Script to do arbutus analysis but just for one body part

library(geiger)
library(arbutus)
library(tidyverse)
library(parallel)

#Load data for brain
path <- "cichlids/expression_data_and_trees/CichlidX_TPM_GeneExpressionMatrix_"
path2 <- "cichlids/expression_data_and_trees/"

br_data <- read.delim(paste0(path, "BR",".txt"))

#Annotate genes vs lnc
dictionary <- read.delim(paste0(path2, "GCF_001858045.1_ASM185804v2_genomic_gtf_gene.txt")) %>%
  select(geneID, biotype)

#Load tree
tree <- read.tree(paste0(path2, "intree"))

#function to separate data into lnc or protein coding
separate_df <- function(matrix, type){
  step1 <- dictionary %>% filter(biotype == type)
  res <- matrix[which(row.names(matrix) %in% step1$geneID),]
  res
}

br_prot <- br_data %>% separate_df("protein_coding") %>% t() %>% treedata(tree, ., sort = TRUE)

#Now split the brain data into pieces
br_df <- br_prot$data 
br_list <- br_df  %>% t() %>% as.data.frame() %>%
  split(rep(1:10, length.out = nrow(df), each = ceiling(nrow(df)/10))) %>%
  map(t)

number <- function( dat_list){
  count = 1
  new_list <- vector(mode = "list", length = length(dat_list))
  for(n in 1:length(dat_list)){
    new_list[[n]] <- dat_list[[n]] %>% list(as.character(n), "prot")
  }
  new_list
}

final_list <- lapply(br_list, function(x) treedata(tree, x, sort = TRUE)) %>% number()

runFC <- function ( td ){
  fitResults <- vector(mode = "list", length = ncol(td$data))
  phy <- td$phy
  data <- td$data
  for(j in 1:ncol(data)){
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

total_process <- function (dat_list){
  avgdat <- dat_list[[1]]
  part <- dat_list[[2]]
  type <- dat_list[[3]]
  fit <- runFC(avgdat)
  fit_name <- paste0("cichlids/arbutus/Fits/brain", part, "_", type)
  saveRDS(fit, file = fit_name)
  df <- model_count(fit)
  aic_name <- paste0("cichlids/arbutus/AIC/AIC_brain_", part, "_", type, ".png")
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  ggsave(aic_name)
  result <- run_arb(fit)
  rds_name <- paste0("cichlids/arbutus/pvals/pvals_brain_", part, "_", type)
  saveRDS(result, file = rds_name)
  result %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("cichlids/arbutus/arbutus_brain_", part, "_", type, ".png")
  ggsave(pval_name)
}

mclapply(final_list, total_process, mc.cores = 10)