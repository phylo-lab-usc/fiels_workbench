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
tree <- read.tree(paste0(path2, "intree")) %>% drop.tip("Calple")

#Make OUwie df
OUwie_reg <- c(rep(1,10), rep(2,7), rep(1,55))
daf <- data.frame(Genus_species = tree$tip.label, Reg = (OUwie_reg))
tree$node.label <- c(rep(1,15), rep(2,6), rep(1,50))
plot(tree, show.node.label = TRUE)

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

runFC <- function ( df){
  fitResults <- vector(mode = "list", length = ncol(df$data))
  phy <- df$phy
  data <- df$data
  phy$node.label <- c(rep(1,15), rep(2,6), rep(1,50))
  for(j in 1:ncol(df$data)){
    fitBM <- tryCatch(fitContinuous(phy, (data[,j]), model = "BM"), error = function(x)list(opt = list(aic = Inf)))
    fitOU <- tryCatch(fitContinuous(phy, (data[,j]), model = "OU"), error = function(x)list(opt = list(aic = Inf)))
    fitEB <- tryCatch(fitContinuous(phy, (data[,j]), model = "EB"), error = function(x)list(opt = list(aic = Inf)))
    OUwie_df <- daf %>% mutate(X = (data[,j]))
    fitBMS <- tryCatch(OUwie(phy, OUwie_df, model = "BMS"), error = function(x)list(AIC = Inf))
    aic <- c(fitBM$opt[["aic"]], fitOU$opt[["aic"]], fitEB$opt[["aic"]], fitBMS$AIC)
    fit <- ifelse(min(aic) == aic[1], list(c(fitBM, model = "BM")), 
                  ifelse(min(aic) == aic[2], list(c(fitOU, model = "OU")), 
                         ifelse(min(aic) == aic[3], list(c(fitEB, model = "EB")),
                                list(c(fitBMS, model = "BMS")))))
    fitResults[j] <- fit
  }
  fitResults
}


model_count <- function (fit) {
  ou = 0
  bm = 0
  eb = 0
  bms = 0
  for(f in fit){
    vec <- f
    ifelse(vec$model == "OU", ou <- ou + 1, ifelse(vec$model == "BM", bm <- bm + 1, ifelse(vec$model == "EB", eb <- eb + 1, bms <- bms + 1)))
  }
  df <- data.frame(OU = ou, BM = bm, EB = eb, BMS = bms)
  b <- df %>% pivot_longer(c(OU, BM, EB, BMS), names_to = "model")
  b
}

class_fix <- function(fitobj){
  if(length(fitobj) == 5) class(fitobj) <- "gfit"
  if(length(fitobj) == 27) class(fitobj) <- "OUwie"
  fitobj
}

run_arb <- function (fits){
  arby <- lapply(fits, function(i) try(arbutus(i), TRUE))
  arby <- arby[sapply(arby, function(x) !inherits(x, "try-error"))]
  arby_df <- map_df(arby, function(i) try(pvalue_arbutus(i), TRUE))
  arby_df
}


total_process <- function (dat_list){
  avgdat <- dat_list[[1]]
  part <- dat_list[[2]]
  type <- dat_list[[3]]
  fit <- runFC(avgdat)
  fit_name <- paste0("cichlids/arbutus/multirate/Fits/brain", part, "_", type)
  saveRDS(fit, file = fit_name)
  df <- model_count(fit)
  aic_name <- paste0("cichlids/arbutus/multirate/AIC/AIC_brain_", part, "_", type, ".png")
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  ggsave(aic_name)
  result <- lapply(fit, class_fix) %>% run_arb()
  rds_name <- paste0("cichlids/arbutus/multirate/pvals/pvals_brain_", part, "_", type)
  saveRDS(result, file = rds_name)
  result %>% select(!m.sig) %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("cichlids/arbutus/multirate/figures/arbutus_brain_", part, "_", type, ".png")
  ggsave(pval_name)
}

#lapply(final_list, total_process)
mclapply(final_list, total_process, mc.cores = 10)

brain1 <- readRDS("cichlids/arbutus/multirate/pvals/pvals_brain_1_prot")
