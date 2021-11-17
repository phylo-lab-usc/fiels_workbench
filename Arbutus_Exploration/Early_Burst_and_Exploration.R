#Comparing and plotting Arbutus exploration with Early Bust model exploration

library(tidyverse)
library(geiger)
library(OUwie)


#Dataset loading for violin plots
EB_violin <- readRDS("Arbutus_Exploration/RDSfiles/EB_data")
others_violin <- readRDS("Arbutus_Exploration/RDSfiles/Exploration3_data")

#Let's just use alpha = 1 for early burst
EB_to_join <- EB_violin %>%
  filter(alpha == "one") %>%
  mutate(model = "EB") %>%
  select(-alpha)

plot_df <- EB_to_join %>% full_join(others_violin)

#Composite violin plot
plot_df %>% pivot_longer(cols = c(-model), names_to = "test.stat") %>%
  ggplot(aes(y = value, x = model, fill = (model))) + 
  geom_violin() + geom_boxplot(width = 0.5) + 
  facet_wrap(~test.stat) + theme_bw()

#Now to show Early Burst fit
EB_EB <- readRDS("Arbutus_Exploration/RDSfiles/EB_fit_to_EB")

#Change value to not a list.
a <- EB_EB %>% select(-c(method, k)) %>% 
  pivot_longer(cols = c(-alpha), names_to = "param") #%>%
  ggplot(aes(x = value)) + geom_boxplot() + theme_classic()

