#Debugging s.hgt script

library(arbutus)
library(geiger)
library(tidyverse)

#Read the data
pvals <- readRDS("Mammal_organs/species_phylogeny/arbutus/pvals_br")
fit <- readRDS("Mammal_organs/species_phylogeny/arbutus/fit_br")

#change fit data to gfit class
gfitClass <- function ( fitObj ) {
  class(fitObj) <- "gfit"
  fitObj
}

retrieve_alpha <- function ( gfit ) {
  model_info(gfit)$pars$alpha
}

retrieve_type <- function ( gfit ) {
  model_info(gfit)$type
}

retrieve_sig <- function ( gfit ) {
  model_info(gfit)$pars$sigsq
}

retrieve_z <- function ( gfit ) {
  model_info(gfit)$pars$z0
}

fit <- fit %>% map(gfitClass)
alpha <- fit %>% map(retrieve_alpha)
type <- fit %>% map(retrieve_type)
sig <- fit %>% map(retrieve_sig)
z <- fit %>% map(retrieve_z)

#4th line of pvals has NA in shgt, compare to 3rd line
#fit1 <- fit[[1]]
#class(fit1) <- "gfit"


#Step by step arbutus
#fits 4-9 are all NA, 1-3 have values
#unit.tree1 <- make_unit_tree(fit1)
#plot(unit.tree1$phy)

#Unit trees for NAs are completely conal, unit trees for others have dimension
#alpha values for all NA trees are around 2.7! Others have alphas that are different.
#Write a script that finds alpha values for all NA values, plots them vs values for non NAs
pvals <- pvals %>% slice_head(n = 100)
test_df <- pvals %>% select(s.hgt) %>% mutate(alpha = alpha, type = type, sig = sig, z = z)
#analysis

#test for if type is related to NA
type_df <- test_df %>% filter(is.na(s.hgt)) #looks like all NAs are OU

#test for NA alpha values
a_na <- test_df %>% filter(is.na(s.hgt)) %>% mutate(alpha = unlist(alpha), group = "NA", sig = unlist(sig), z = unlist(z))
alpha_not_na <- test_df %>% filter(!is.na(s.hgt)) %>% filter(type == "OU") %>% mutate(alpha = unlist(alpha), group = "Not", sig = unlist(sig), z = unlist(z))

all_df <- full_join(a_na, alpha_not_na)
all_df %>% ggplot(aes(y = alpha, x = group, color = group)) + geom_point(position = position_jitter(width = 0.1)) + theme_bw()
all_df %>% ggplot(aes(y = sig, x = group, color = group)) + geom_point(position = position_jitter(width = 0.1))
all_df %>% ggplot(aes(y = z, x = group, color = group)) + geom_point(position = position_jitter(width = 0.1))

#It looks like a BOUNDS error
unit_list <- map(fit, make_unit_tree)
obs <- calculate_pic_stat(unit_list)
sim.dat <- simulate_char_unit(unit_list)
sim <- calculate_pic_stat(sim.dat)
res <- compare_pic_stat(obs, sim)
plot(res)
pvalue_arbutus(res)