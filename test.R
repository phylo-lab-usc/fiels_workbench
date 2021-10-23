#Test script

library(geiger)
library(arbutus)
library(ape)

data(finch)
geo <- finch$phy
dat <- finch$data[,"wingL"]

first_sim <- sim.char(geo, 0.02, 1000)

unit_tree <- make_unit_tree(geo, data = dat)

first_sim_unit <- simulate_char_unit(unit_tree)

test.stat <- calculate_pic_stat(unit_tree)
