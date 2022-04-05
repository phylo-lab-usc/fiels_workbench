#temp script

#Arbutus analysis

library(geiger)
library(arbutus)
library(tidyverse)
library(parallel)

#Load each organ
path <- "cichlids/expression_data_and_trees/CichlidX_TPM_GeneExpressionMatrix_"
path2 <- "cichlids/expression_data_and_trees/"

br_data <- read.delim(paste0(path, "BR",".txt"))
lp_data <- read.delim(paste0(path, "LP",".txt"))
te_data <- read.delim(paste0(path, "TE",".txt"))
ov_data <- read.delim(paste0(path, "OV",".txt"))
ve_data <- read.delim(paste0(path, "VE",".txt"))
gi_data <- read.delim(paste0(path, "GI",".txt"))

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

br_lnc <- br_data %>% separate_df("lncRNA") %>% t() %>% treedata(tree, ., sort = TRUE)
br_prot <- br_data %>% separate_df("protein_coding") %>% t() %>% treedata(tree, ., sort = TRUE)

lp_lnc <- lp_data %>% separate_df("lncRNA") %>% t() %>% treedata(tree, ., sort = TRUE)
lp_prot <- lp_data %>% separate_df("protein_coding") %>% t() %>% treedata(tree, ., sort = TRUE)

te_lnc <- te_data %>% separate_df("lncRNA") %>% t() %>% treedata(tree, ., sort = TRUE)
te_prot <- te_data %>% separate_df("protein_coding") %>% t() %>% treedata(tree, ., sort = TRUE)

ov_lnc <- ov_data %>% separate_df("lncRNA") %>% t() %>% treedata(tree, ., sort = TRUE)
ov_prot <- ov_data %>% separate_df("protein_coding") %>% t() %>% treedata(tree, ., sort = TRUE)

ve_lnc <- ve_data %>% separate_df("lncRNA") %>% t() %>% treedata(tree, ., sort = TRUE)
ve_prot <- ve_data %>% separate_df("protein_coding") %>% t() %>% treedata(tree, ., sort = TRUE)

gi_lnc <- gi_data %>% separate_df("lncRNA") %>% t() %>% treedata(tree, ., sort = TRUE)
gi_prot <- gi_data %>% separate_df("protein_coding") %>% t() %>% treedata(tree, ., sort = TRUE)

rm(br_data, lp_data, te_data, ov_data, ve_data, gi_data, dictionary)

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
  fit_name <- paste0("cichlids/arbutus/Fits/", part, "_", type)
  saveRDS(fit, file = fit_name)
  df <- model_count(fit)
  aic_name <- paste0("cichlids/arbutus/AIC/AIC_", part, "_", type, ".png")
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  ggsave(aic_name)
  result <- run_arb(fit)
  rds_name <- paste0("cichlids/arbutus/pvals/pvals_", part, "_", type)
  saveRDS(result, file = rds_name)
  result %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("cichlids/arbutus/arbutus_", part, "_", type, ".png")
  ggsave(pval_name)
}

br_list1 <- list(br_lnc, "br", "lnc")
br_list2 <- list(br_prot, "br", "prot")
lp_list1 <- list(lp_lnc, "lp", "lnc")
lp_list2 <- list(lp_prot, "lp", "prot")
te_list1 <- list(te_lnc, "te", "lnc")
te_list2 <- list(te_prot, "te", "prot")
ov_list1 <- list(ov_lnc, "ov", "lnc")
ov_list2 <- list(ov_prot, "ov", "prot")
ve_list1 <- list(ve_lnc, "ve", "lnc")
ve_list2 <- list(ve_prot, "ve", "prot")
gi_list1 <- list(gi_lnc, "gi", "lnc")
gi_list2 <- list(gi_prot, "gi", "prot")

all_list <- list(br_list2, lp_list2,
                 te_list2, ov_list2,
                 ve_list2, gi_list2)

mclapply(all_list, total_process, mc.cores = 6)
