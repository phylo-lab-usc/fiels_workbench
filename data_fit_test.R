#Testing fitContinuous with OU models

#First, load required packages

library(tidyverse)
library(ape)
library(geiger)
library(OUwie)
library(broom)

#Now, simulate birth-death tree

tr <- sim.bdtree(n = 128)
plot(tr)

#First test: Simulate character data once, and then fitContinuous 1000 times, with and without fixing alpha bounds

#simulate data use rTraitCont
char_df <- data.frame(rTraitCont(tr, model = "OU"))

#not fixing bounds
OU_fit1 <- replicate(10, fitContinuous(tr, char_df, model = "OU"))  

#fixing bounds
OU_fit2 <- replicate(10, fitContinuous(tr, char_df, model = "OU", bounds = list(alpha = c(0,1.1))))

#What we learned: fitContinuous gets the same result every time, so need to simulate character data

#Second test: Simulate character data 1000 times, then fit continuous to all
sim_and_fit <- function (tree) {
  df <- data.frame(rTraitCont(tree, model = "OU"))
  fit <- fitContinuous(tr, df, model = "OU")
  fit$opt
}

char_list <- replicate(1000, sim_and_fit(tr))
char_df2 <- fix_data_frame(t(data.frame(char_list)))
char_df2$alpha <- unlist(char_df2$alpha)
char_df2 %>% ggplot(aes(x = alpha)) + geom_boxplot() + theme_classic() + annotate("text", x = 2, y = -0.2, label = mean(char_df2$alpha))
