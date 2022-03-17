#temp script to run arbutus

library(tidyverse)
library(parallel)
library(arbutus)
library(geiger)

str_list <- list.files("Mammal_organs/gene_family_phylogenies/arbutus/fits/")

run_arb <- function (fits){
  arby <- vector("list", length = length(fits))
  count = 1
  for(f in fits){
    class(f) <- "gfit"
    arby[[count]] <- try(arbutus(f))
    count = count + 1
  }
  arby <- arby %>% purrr::discard(~ !typeof(.x) == "list")
  arby_df <- map_df(arby, pvalue_arbutus)
  arby_df
}

total_process <- function (string){
  fit <- readRDS(paste0("Mammal_organs/gene_family_phylogenies/arbutus/fits/", string))
  part <- str_sub(string, 5,6)
  result <- run_arb(fit)
  rds_name <- paste0("Mammal_organs/gene_family_phylogenies/arbutus/pvals/pvals_", part)
  saveRDS(result, file = rds_name)
  result %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("Mammal_organs/gene_family_phylogenies/arbutus/figures/arbutus_", part, ".png")
  ggsave(pval_name)
}


#total_process("fit_br")
mclapply(str_list, total_process, mc.cores = 6)
