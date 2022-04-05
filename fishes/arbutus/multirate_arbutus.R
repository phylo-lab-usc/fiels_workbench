#Analysis including multi-rate BMS

library(geiger)
library(arbutus)
library(tidyverse)
library(parallel)

#Load the data
tree <- read.tree("fishes/data/recodedTreeNamed.tre") %>% chronoMPL()

OUwie_reg <- c(1,2,2,1,2,1,2,1,2,1,1,2,1,2,2,1,1,2,1,2)
df <- data.frame(Genus_species = tree$tip.label, Reg = OUwie_reg)

data_ave <- read.csv("fishes/data/master_fpkm.csv", row.names = 1) %>% rownames_to_column("genes") %>% group_by(genes)%>%
  transmute(GholNS = mean(c(GholNS_1, GholNS_2, GholNS_3, GholNS_4, GholNS_5, GholNS_6)),
            GsexNS = mean(c(GsexNS_1, GsexNS_2, GsexNS_3, GsexNS_4, GsexNS_5, GsexNS_6)),
            PbimNS = mean(c(PbimNS_1, PbimNS_2, PbimNS_3, PbimNS_4, PbimNS_5, PbimNS_6)),
            XhelNS = mean(c(XhelNS_1, XhelNS_2, XhelNS_3, XhelNS_4, XhelNS_5, XhelNS_6)),
            LperNS = mean(c(LperNS_1, LperNS_2, LperNS_3, LperNS_4, LperNS_5, LperNS_6)),
            PlatNS = mean(c(PlatNS_1, PlatNS_2, PlatNS_3, PlatNS_4, PlatNS_5, PlatNS_6)),
            PlimNS = mean(c(PlimNS_1, PlimNS_2, PlimNS_3, PlimNS_4, PlimNS_5, PlimNS_6)),
            PmexPuyNS = mean(c(PmexPuyNS_1, PmexPuyNS_2, PmexPuyNS_3, PmexPuyNS_4, PmexPuyNS_5, PmexPuyNS_6)),
            PmexPichNS = mean(c(PmexPichNS_1, PmexPichNS_2, PmexPichNS_3, PmexPichNS_4, PmexPichNS_5, PmexPichNS_6)),
            PmexTacNS = mean(c(PmexTacNS_1, PmexTacNS_2, PmexTacNS_3, PmexTacNS_4, PmexTacNS_5, PmexTacNS_6)),
            PmexPuyS = mean(c(PmexPuyS_1, PmexPuyS_2, PmexPuyS_3, PmexPuyS_4, PmexPuyS_5)),
            PmexTacS = mean(c(PmexTacS_2, PmexTacS_3, PmexTacS_4, PmexTacS_5, PmexTacS_6)),
            PmexPichS = mean(c(PmexPichS_1, PmexPichS_2, PmexPichS_3, PmexPichS_4, PmexPichS_5, PmexPichS_6)),
            PlatS = mean(c(PlatS_1, PlatS_2, PlatS_3, PlatS_4, PlatS_5, PlatS_6)),
            LsulS = mean(c(LsulS_1, LsulS_2, LsulS_3, LsulS_4, LsulS_5, LsulS_6)),
            GsexS = mean(c(GsexS_1, GsexS_2, GsexS_3, GsexS_4, GsexS_5, GsexS_6)),
            GeurS = mean(c(GeurS_1, GeurS_2, GeurS_3, GeurS_4, GeurS_5, GeurS_6)),
            GholS = mean(c(GholS_1, GholS_2, GholS_3, GholS_4, GholS_5, GholS_6)),
            PbimS = mean(c(PbimS_1, PbimS_2, PbimS_3, PbimS_4, PbimS_5, PbimS_6)),
            XhelS = mean(c(XhelS_1, XhelS_2, XhelS_3, XhelS_4, XhelS_5, XhelS_6))) %>% column_to_rownames("genes") %>%
  split(rep(1:10, length.out = nrow(df), each = ceiling(nrow(df)/10))) %>%
  map(t)

standard_error <- function(x) sd(x) / sqrt(length(x))

data_SE <- read.csv("fishes/data/master_fpkm.csv", row.names = 1) %>% rownames_to_column("genes") %>% group_by(genes) %>%
  summarise(genes, SE = standard_error(across(GholNS_1:XhelS_6))) %>% 
  split(rep(1:10, length.out = nrow(df), each = ceiling(nrow(df)/10)))

#Functions
number <- function( dat_list, SE_list ){
  count = 1
  new_list <- vector(mode = "list", length = length(dat_list))
  for(n in 1:length(dat_list)){
    new_list[[n]] <- dat_list[[n]] %>% list(as.character(n), SE_list[[n]])
  }
  new_list
}

runFC <- function ( df, StE ){
  fitResults <- vector(mode = "list", length = ncol(df))
  td <- treedata(tree, df, sort = TRUE)
  phy <- td$phy
  data <- td$data
  phy$node.label <- c(1,1,1,1,2,2,2,2,2,1,1,1,2,2,1,1,2,2,2)
  for(j in 1:ncol(df)){
    fitBM <- fitContinuous(phy, data[,j], StE[[2]][[j]], model = "BM")
    fitOU <- fitContinuous(phy, data[,j], StE[[2]][[j]], model = "OU")
    fitEB <- fitContinuous(phy, data[,j], StE[[2]][[j]], model = "EB")
    OUwie_df <- df %>% mutate(X = data[,j])
    fitBMS <- tryCatch(OUwie(phy, OUwie_df, model = "BMS"), error = function(x)NULL)
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

run_arb <- function (fits){
  arby <- vector("list", length = length(fits))
  count = 1
  for(f in fits){
    arby[[count]] <- arbutus(f)
    count = count + 1
  }
  arby_df <- map_df(arby, pvalue_arbutus)
  arby_df
}

total_process <- function (dat_obj){
  dat_matrix <- dat_obj[[1]]
  num <- dat_obj[[2]]
  se <- dat_obj[[3]]
  fit <- runFC(dat_matrix, se)
  fit_name <- paste0("fishes/arbutus/multi-rate/Fits/fit_", num)
  saveRDS(fit, file = fit_name)
  df <- model_count(fit)
  aic_name <- paste0("fishes/arbutus/multi-rate/AIC/AIC_", num, ".png")
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  ggsave(aic_name)
  result <- run_arb(fit)
  rds_name <- paste0("fishes/arbutus/multi-rate/pvals/pvals_", num)
  saveRDS(result, file = rds_name)
  result %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("fishes/arbutus/multi-rate/plots/arbutus_", num, ".png")
  ggsave(pval_name)
}


data_modified <- number(data_ave, data_SE)
mclapply(data_modified, total_process, mc.cores = 10)