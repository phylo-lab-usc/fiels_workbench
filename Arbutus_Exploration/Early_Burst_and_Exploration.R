#Comparing and plotting Arbutus exploration with Early Bust model exploration

library(tidyverse)
library(geiger)
library(OUwie)



#Now to show Early Burst fit
EB_EB <- readRDS("Arbutus_Exploration/RDSfiles/EB_fit_to_EB")
EB_BM <- readRDS("Arbutus_Exploration/RDSfiles/EB_fit_to_BM")
EB_OU <- readRDS("Arbutus_Exploration/RDSfiles/EB_fit_to_OU")

#Plot fit of different models to rescaled EB tree.
EB_EB %>% select(-c(method, k)) %>% 
  pivot_longer(cols = c(-alpha), names_to = "param") %>%
  mutate(value = unlist(value)) %>%
  ggplot(aes(x = value, fill = alpha)) +
  geom_boxplot() + facet_wrap(~param, scales = "free") +
  theme_classic() + labs(title = "Parameters for Early Burst Model fit to Scaled Tree")

EB_BM %>% select(-c(method, k)) %>% 
  pivot_longer(cols = c(-alpha), names_to = "param") %>%
  mutate(value = unlist(value)) %>%
  ggplot(aes(x = value, fill = alpha)) + geom_boxplot() + facet_wrap(~param, scales = "free") + theme_classic() 

EB_OU %>% select(-c(method, k)) %>% 
  pivot_longer(cols = c(-alpha), names_to = "param") %>%
  mutate(value = unlist(value)) %>%
  ggplot(aes(x = value, fill = alpha)) + geom_boxplot() + facet_wrap(~param, scales = "free") + theme_classic() 
