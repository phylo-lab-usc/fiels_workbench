library(biomaRt)
library(tidyverse)
library(parallel)

runlist <- list.files("Mammal_organs/gene_family_phylogenies/fasta_files/") %>% map(sub3) %>% as.character()
seqlist <- ortho %>% filter(! Human %in% runlist) %>% t() %>% as.data.frame() %>% as.list()

mclapply(seqlist, get_sequences_list, mc.cores = 12)
