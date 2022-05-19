#Analysis of snake venom toxins

library(geiger)
library(arbutus)
library(tidyverse)
library(parallel)

#Load the data
tree <- read.nexus("snakes/data/Final_tree.nex")
data <- read.csv("snakes/data/vPhen_data.csv") %>% column_to_rownames(var = "species") %>%
  select(-c(Total, cm, Family, Author, Tech, Ref)) %>% as.matrix()

#Make OUwie df
#OUwie_reg <- c(rep(1,10), rep(2,7), rep(1,55))
#daf <- data.frame(Genus_species = tree$tip.label, Reg = (OUwie_reg))
#tree$node.label <- c(rep(1,15), rep(2,6), rep(1,50))

#Functions
runFC <- function ( df ){
  fitResults <- vector(mode = "list", length = ncol(df))
  for(j in 1:ncol(df)){
    td <- treedata(tree, data[,j], sort = TRUE)
    #td$phy$node.label <- 
    fitBM <- fitContinuous(td$phy, td$data, model = "BM")
    fitOU <- fitContinuous(td$phy, td$data, model = "OU")
    fitEB <- fitContinuous(td$phy, td$data, model = "EB")
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

total_process <- function (dat_matrix){
  fit <- runFC(dat_matrix)
  fit_name <- paste0("snakes/arbutus/BMS/fit")
  saveRDS(fit, file = fit_name)
  df <- model_count(fit)
  aic_name <- paste0("snakes/arbutus/BMS/AIC.png")
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  ggsave(aic_name)
  result <- lapply(fit, class_fix) %>% run_arb()
  rds_name <- paste0("snakes/arbutus/BMS/pvals")
  saveRDS(result, file = rds_name)
  result %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("snakes/arbutus/BMS/arbutus.png")
  ggsave(pval_name)
}

total_process(data)