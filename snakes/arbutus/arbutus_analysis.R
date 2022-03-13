#Analysis of snake venom toxins

library(geiger)
library(arbutus)
library(tidyverse)

#Load the data
tree <- read.nexus("snakes/data/Final_tree.nex")
data <- read.csv("snakes/data/vPhen_data.csv") %>% column_to_rownames(var = "species") %>%
  select(-c(Total, cm, Family, Author, Tech, Ref)) %>% as.matrix()

#Functions
runFC <- function ( df ){
  fitResults <- vector(mode = "list", length = ncol(df))
  for(j in 1:ncol(df)){
    td <- treedata(tree, data[,j], sort = TRUE)
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

total_process <- function (dat_matrix){
  fit <- runFC(dat_matrix)
  fit_name <- paste0("snakes/arbutus/fit")
  saveRDS(fit, file = fit_name)
  df <- model_count(fit)
  aic_name <- paste0("snakes/arbutus/AIC.png")
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  ggsave(aic_name)
  result <- run_arb(fit)
  rds_name <- paste0("snakes/arbutus/pvals")
  saveRDS(result, file = rds_name)
  result %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("snakes/arbutus/arbutus.png")
  ggsave(pval_name)
}

total_process(data)