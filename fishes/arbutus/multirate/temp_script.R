#temp_script

library(geiger)
library(arbutus)
library(tidyverse)
library(parallel)
library(OUwie)

path = "fishes/arbutus/multirate/"
fit_names <- list.files(paste0(path, "Fits"))

readFits <- function(name){
  res <- readRDS(paste0(path, "Fits/", name))
  list(res, substring(name,5))
}

fit_list <- lapply(fit_names, readFits)

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

temp_process <- function (list_obj){
  fit <- list_obj[[1]]
  num <- list_obj[[2]]
  df <- model_count(fit)
  aic_name <- paste0("fishes/arbutus/multirate/AIC/AIC_", num, ".png")
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  ggsave(aic_name)
  result <- run_arb(fit)
  rds_name <- paste0("fishes/arbutus/multirate/pvals/pvals_", num)
  saveRDS(result, file = rds_name)
  result %>% select(!m.sig) %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("fishes/arbutus/multirate/plots/arbutus_", num, ".png")
  ggsave(pval_name)
}


mclapply(fit_list, temp_process, mc.cores = 20)