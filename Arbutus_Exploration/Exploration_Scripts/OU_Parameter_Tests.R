#OU parameter testing

#First, load required packages

library(tidyverse)
library(ape)
library(geiger)
library(OUwie)
library(broom)
library(arbutus)

#Now, simulate birth-death tree

tr <- sim.bdtree(n = 128)

dat1 <- data.frame(tr$tip.label) %>% mutate(Reg = c(rep(1, 64), rep(2, 64))) %>% rename(Genus_species = tr.tip.label)
dat2 <- data.frame(tr$tip.label) %>% mutate(Reg = 1) %>% rename(Genus_species = tr.tip.label)

#First test parameters for regular OU
tree_OU <- tr
tree_OU$node.label <- rep(1,127)

sim_and_fit_OU <- function (tree, dat, alpha, sigma, theta) {
  df <- OUwie.sim(tree, data = dat, alpha = alpha, sigma.sq = sigma, theta0 = 1, theta = theta)
  df_fix <- df
  row.names(df_fix) <- df_fix$Genus_species
  df_fix <- df_fix %>% select(X)
  fit <- fitContinuous(tree, df_fix, model = "OU")
  a <- arbutus(fit)
  a
}

#Test different alpha. Alpha term regulates "strength" of optima
#Single alpha: Small, medium, large
alpha.sm <- replicate(1000, sim_and_fit_OU(tree_OU, dat2, c(0.01, 0.01), c(0.9,0.9), c(1.0, 1.0)))
alpha.med <- replicate(1000, sim_and_fit_OU(tree_OU, dat2, c(0.1, 0.1), c(0.9,0.9), c(1.0, 1.0)))
alpha.lar <- replicate(1000, sim_and_fit_OU(tree_OU, dat2, c(1.0, 1.0), c(0.9,0.9), c(1.0, 1.0)))

#Test different sigma. Sigma.sq tests evolutionary rate. 
#Single rate
sigma.sm <- replicate(1000, sim_and_fit_OU(tree_OU, dat2, c(0.1, 0.1), c(0.1,0.1), c(1.0, 1.0)))
sigma.med <- replicate(1000, sim_and_fit_OU(tree_OU, dat2, c(0.1, 0.1), c(0.5,0.5), c(1.0, 1.0)))
sigma.lar <- replicate(1000, sim_and_fit_OU(tree_OU, dat2, c(0.1, 0.1), c(1.0,1.0), c(1.0, 1.0)))

#Test different theta. Theta represents evolutionary optima. Root theta is 1
theta.sm <- replicate(1000, sim_and_fit_OU(tree_OU, dat2, c(0.1, 0.1), c(0.9,0.9), c(1.1, 1.1)))
theta.med <- replicate(1000, sim_and_fit_OU(tree_OU, dat2, c(0.1, 0.1), c(0.9,0.9), c(1.5, 1.5)))
theta.lar <- replicate(1000, sim_and_fit_OU(tree_OU, dat2, c(0.1, 0.1), c(0.9,0.9), c(2.0, 2.0)))

#Now test for multivariate OU
tree_OUwie <- tr
tree_OUwie$node.label <- c(rep(1, 64), rep(2, 63))

#Multiple alpha: Small diff, medium diff, large diff
alpha.smdiff <- replicate(1000, sim_and_fit_OU(tree_OUwie, dat1, c(1.0, 0.95), c(0.9,0.9), c(1.0, 2.0)))
alpha.meddiff <- replicate(1000, sim_and_fit_OU(tree_OUwie, dat1, c(1.0, 0.5), c(0.9,0.9), c(1.0, 2.0)))
alpha.lardiff <- replicate(1000, sim_and_fit_OU(tree_OUwie, dat1, c(1.0, 0.05), c(0.9,0.9), c(1.0, 2.0)))

#Multisigma
sigma.smdiff <- replicate(1000, sim_and_fit_OU(tree_OUwie, dat1, c(0.1, 0.1), c(0.45,0.5), c(1.0, 2.0)))
sigma.meddiff <- replicate(1000, sim_and_fit_OU(tree_OUwie, dat1, c(0.1, 0.1), c(0.9,0.45), c(1.0, 2.0)))
sigma.lardiff <- replicate(1000, sim_and_fit_OU(tree_OUwie, dat1, c(0.1, 0.1), c(0.9,0.1), c(1.0, 2.0)))

#Multitheta
theta.smdiff <- replicate(1000, sim_and_fit_OU(tree_OUwie, dat1, c(0.1, 0.1), c(0.45,0.45), c(1.0, 1.5)))
theta.meddiff <- replicate(1000, sim_and_fit_OU(tree_OUwie, dat1, c(0.1, 0.1), c(0.45,0.45), c(1.0, 2.0)))
theta.lardiff <- replicate(1000, sim_and_fit_OU(tree_OUwie, dat1, c(0.1, 0.1), c(0.45,0.45), c(1.0, 3.0)))

#retrieve pvals
#alpha
alpha.sm_pvals <- alpha.sm[1,]
alpha.med_pvals <- alpha.med[1,]
alpha.lar_pvals <- alpha.lar[1,]
alpha.smdiff_pvals <- alpha.smdiff[1,]
alpha.meddiff_pvals <- alpha.meddiff[1,]
alpha.lardiff_pvals <- alpha.lardiff[1,]

#sigma
sigma.sm_pvals <- sigma.sm[1,]
sigma.med_pvals <- sigma.med[1,]
sigma.lar_pvals <- sigma.lar[1,]
sigma.smdiff_pvals <- sigma.smdiff[1,]
sigma.meddiff_pvals <- sigma.meddiff[1,]
sigma.lardiff_pvals <- sigma.lardiff[1,]

#theta
theta.sm_pvals <- theta.sm[1,]
theta.med_pvals <- theta.med[1,]
theta.lar_pvals <- theta.lar[1,]
theta.smdiff_pvals <- theta.smdiff[1,]
theta.meddiff_pvals <- theta.meddiff[1,]
theta.lardiff_pvals <- theta.lardiff[1,]

#arbutus_transform() from custom_functions.R

#transform
alpha.sm_df <- arbutus_transform(alpha.sm_pvals, 1000) %>% mutate(size = "sm")
alpha.med_df <- arbutus_transform(alpha.med_pvals, 1000) %>% mutate(size = "med")
alpha.lar_df <- arbutus_transform(alpha.lar_pvals, 1000) %>% mutate(size = "lar")
alpha.smdiff_df <- arbutus_transform(alpha.smdiff_pvals, 1000) %>% mutate(size = "smdiff")
alpha.meddiff_df <- arbutus_transform(alpha.meddiff_pvals, 1000) %>% mutate(size = "meddiff")
alpha.lardiff_df <- arbutus_transform(alpha.lardiff_pvals, 1000) %>% mutate(size = "lardiff")

sigma.sm_df <- arbutus_transform(sigma.sm_pvals, 1000) %>% mutate(size = "sm")
sigma.med_df <- arbutus_transform(sigma.med_pvals, 1000) %>% mutate(size = "med")
sigma.lar_df <- arbutus_transform(sigma.lar_pvals, 1000) %>% mutate(size = "lar")
sigma.smdiff_df <- arbutus_transform(sigma.smdiff_pvals, 1000) %>% mutate(size = "smdiff")
sigma.meddiff_df <- arbutus_transform(sigma.meddiff_pvals, 1000) %>% mutate(size = "meddiff")
sigma.lardiff_df <- arbutus_transform(sigma.lardiff_pvals, 1000) %>% mutate(size = "lardiff")

theta.sm_df <- arbutus_transform(theta.sm_pvals, 1000) %>% mutate(size = "sm")
theta.med_df <- arbutus_transform(theta.med_pvals, 1000) %>% mutate(size = "med")
theta.lar_df <- arbutus_transform(theta.lar_pvals, 1000) %>% mutate(size = "lar")
theta.smdiff_df <- arbutus_transform(theta.smdiff_pvals, 1000) %>% mutate(size = "smdiff")
theta.meddiff_df <- arbutus_transform(theta.meddiff_pvals, 1000) %>% mutate(size = "meddiff")
theta.lardiff_df <- arbutus_transform(theta.lardiff_pvals, 1000) %>% mutate(size = "lardiff")

#fuse
alpha_df <- full_join(alpha.sm_df, alpha.med_df) %>% full_join(alpha.lar_df) %>%
  full_join(alpha.smdiff_df) %>% full_join(alpha.meddiff_df) %>% full_join(alpha.lardiff_df)
sigma_df <- full_join(sigma.sm_df, sigma.med_df) %>% full_join(sigma.lar_df) %>%
  full_join(sigma.smdiff_df) %>% full_join(sigma.meddiff_df) %>% full_join(sigma.lardiff_df)
theta_df <- full_join(theta.sm_df, theta.med_df) %>% full_join(theta.lar_df) %>%
  full_join(theta.smdiff_df) %>% full_join(theta.meddiff_df) %>% full_join(theta.lardiff_df)

#save data
saveRDS(alpha_df, "Arbutus_Exploration/RDSfiles/OU_alpha_data")
saveRDS(sigma_df, "Arbutus_Exploration/RDSfiles/OU_sigma_data")
saveRDS(theta_df, "Arbutus_Exploration/RDSfiles/OU_theta_data")


#pivot and plot
alpha_df %>% pivot_longer(cols = c(-size), names_to = "test.stat") %>%
  ggplot(aes(y = value, x = size, fill = (size))) + geom_violin() + geom_boxplot(width = 0.5) + facet_wrap(~test.stat) + theme_bw()
ggsave("Arbutus_Exploration/OU_alpha_plot.png")

sigma_df %>% pivot_longer(cols = c(-size), names_to = "test.stat") %>%
  ggplot(aes(y = value, x = size, fill = (size))) + geom_violin() + geom_boxplot(width = 0.5) + facet_wrap(~test.stat) + theme_bw()
ggsave("Arbutus_Exploration/OU_sigma_plot.png")

theta_df %>% pivot_longer(cols = c(-size), names_to = "test.stat") %>%
  ggplot(aes(y = value, x = size, fill = (size))) + geom_violin() + geom_boxplot(width = 0.5) + facet_wrap(~test.stat) + theme_bw()
ggsave("Arbutus_Exploration/OU_theta_plot.png")
