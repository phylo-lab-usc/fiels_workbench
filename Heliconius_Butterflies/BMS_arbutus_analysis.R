#Multirate Arbutus Analysis


library(geiger)
library(arbutus)
library(tidyverse)
library(flipR)
library(parallel)

#Load the data
tr <- read.nexus("Heliconius_Butterflies/Data/Heliconiini.multiple_uniform_constraints.MCC.tree")
plot(tr)

tree <- keep.tip(tr, c("Heliconius_charithonia_8830_P", "Heliconius_sara_8862_P",
                       "Heliconius_erato_erato_NCS2556_FG", "Heliconius_doris_02_1939_Pe",
                       "Heliconius_melpomene_melpomene_9317_FG"))
#Rename tips
tree$tip.label <- c("Heliconius_charithonia", "Heliconius_doris", "Heliconius_erato",
                    "Heliconius_melpomene", "Heliconius_sara")

plot(tree)

exp_matr <- read.csv("Heliconius_Butterflies/Data/expression_matrix.csv") #%>%

#Split into sexes
female_dat <- exp_matr %>% as.data.frame() %>% select(Orthogroups, contains("_F"))
male_dat <- exp_matr%>% as.data.frame() %>% select(Orthogroups, contains("_M"))

#Take averages for each species by sex
female_avg_dat <- female_dat %>% group_by(Orthogroups) %>%
  transmute("Heliconius_charithonia" = mean(HCH_1_F, HCH_2_F, HCH_3_F, HCH_4_F, HCH_5_F, HCH_6_F),
            "Heliconius_doris" = mean(HDO_1_F, HDO_2_F, HDO_3_F, HDO_4_F, HDO_5_F, HDO_6_F),
            "Heliconius_erato" = mean(HER_1_F, HER_2_F, HER_3_F),
            "Heliconius_melpomene" = mean(HME_1_F, HME_2_F, HME_3_F, HME_4_F),
            "Heliconius_sara" = mean(HSA_1_F, HSA_2_F, HSA_3_F, HSA_4_F, HSA_5_F))

male_avg_dat <- male_dat %>% group_by(Orthogroups) %>%
  transmute("Heliconius_charithonia" = mean(HCH_1_M, HCH_2_M, HCH_3_M, HCH_4_M, HCH_5_M, HCH_6_M),
            "Heliconius_doris" = mean(HDO_1_M, HDO_2_M, HDO_3_M, HDO_4_M, HDO_5_M, HDO_6_M),
            "Heliconius_erato" = mean(HER_1_M, HER_2_M, HER_3_M),
            "Heliconius_melpomene" = mean(HME_1_M, HME_2_M, HME_3_M, HME_4_M),
            "Heliconius_sara" = mean(HSA_1_M, HSA_2_M, HSA_3_M, HSA_4_M, HSA_5_M))

#Standard Error function
standard_error <- function(x) sd(x) / sqrt(length(x))

#Get SE for males and females
female_SE <- female_dat %>% group_by(Orthogroups) %>%
  summarise(Orthogroups, SE = standard_error(across(HCH_1_F:HSA_5_F)))

male_SE <- male_dat %>% group_by(Orthogroups) %>%
  summarise(Orthogroups, SE = standard_error(across(HCH_1_M:HSA_5_M)))

#Remove unnecessary data
rm(tr, male_dat, female_dat, exp_matr)

#Now need to flip tables and properly format
format_expr_data <- function (avgdat) {
  temp <- avgdat %>% pull(Orthogroups)
  avgdat <- avgdat %>% ungroup() %>% select(!Orthogroups)
  dat <- flip(avgdat)
  colnames(dat) <- temp
  dat
}

#Running fitcontinuous

runFC <- function ( dat, SE ){
  fitResults <- vector(mode = "list", length = ncol(dat))
  tdf <- treedata(tree, dat, sort = TRUE)
  phy <- tdf$phy
  data <- tdf$data
  for(j in 1:ncol(dat)){
    fitBM <- fitContinuous(phy, data[,j], SE[[2]][[j]], model = "BM")
    fitOU <- fitContinuous(phy, data[,j], SE[[2]][[j]], model = "OU")
    fitEB <- fitContinuous(phy, data[,j], SE[[2]][[j]], model = "EB")
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

total_process <- function (dat_list){
  avgdat <- dat_list[[1]]
  part <- dat_list[[2]]
  SE <- dat_list[[3]]
  exp <- format_expr_data(avgdat)
  fit <- runFC(exp, SE)
  fit_name <- paste0("Heliconius_Butterflies/arbutus/fit_", part)
  saveRDS(fit, file = fit_name)
  df <- model_count(fit)
  aic_name <- paste0("Heliconius_Butterflies/arbutus/AIC_", part, ".png")
  df %>% ggplot(aes(model, value)) + geom_col() + theme_classic()
  ggsave(aic_name)
  result <- run_arb(fit)
  rds_name <- paste0("Heliconius_Butterflies/arbutus/pvals_", part)
  saveRDS(result, file = rds_name)
  result %>% pivot_longer(cols = everything(), names_to = "tstat") %>%
    ggplot(aes(value)) + geom_histogram(aes(y = ..density..)) + facet_wrap(~tstat, nrow = 1) + theme_bw()
  pval_name <- paste0("Heliconius_Butterflies/arbutus/arbutus_", part, ".png")
  ggsave(pval_name)
}

female_list <- list(female_avg_dat, "female", female_SE)
male_list <- list(male_avg_dat, "male", male_SE)

all_list <- list(female_list, male_list)

mclapply(all_list, total_process, mc.cores = 4)