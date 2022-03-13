#Analysis of extremophile fish

library(geiger)
library(arbutus)
library(tidyverse)
library(parallel)

#Load the data
tree <- read.tree("fishes/data/recodedTreeNamed.tre")
data <- read.csv("fishes/data/master_fpkm.csv", row.names = 1) %>%
  split(rep(1:10, length.out = nrow(df), each = ceiling(nrow(df)/10))) %>%
  map(t)

#Functions
number <- function( dat_list ){
  count = 1
  new_list <- vector(mode = "list", length = length(dat_list))
  for(n in 1:length(dat_list)){
    new_list[[n]] <- dat_list[[n]] %>% list(as.character(n))
  }
  new_list
}

runFC <- function ( df ){
  fitResults <- vector(mode = "list", length = ncol(df))
  for(j in 1:ncol(df)){
    td <- treedata(tree, df[,j])
    fitBM <- fitContinuous(td$phy, td$data, model = "BM")
    fitOU <- fitContinuous(td$phy, td$data, model = "OU")
    fitEB <- fitContinuous(td$phy, td$data, model = "EB")
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

total_process <- function (dat_obj){
  dat_matrix <- dat_obj[[1]]
  num <- dat_obj[[2]]
  fit <- runFC(dat_matrix)
  fit_name <- paste0("fishes/arbutus/Fits/fit_", num)
  saveRDS(fit, file = fit_name)
  df <- model_count(fit)
  aic_name <- paste0("fishes/arbutus/AIC/AIC", num, ".png")
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  ggsave(aic_name)
  result <- run_arb(fit)
  rds_name <- paste0("fishes/arbutus/pvals/pvals", num)
  saveRDS(result, file = rds_name)
  result %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("fishes/arbutus/plots/arbutus", num, ".png")
  ggsave(pval_name)
}


data_modified <- number(data)

mclapply(data_modified, total_process, mc.cores = 10)