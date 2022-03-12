#Get sequences for each gene ortholog group via biomart

library(biomaRt)
library(tidyverse)
library(parallel)

ortho <- read.delim("Mammal_organs/Supplementary_Data2/Ortho_1to1_AllSpecies.txt", sep = " ")
seqlist <- ortho %>% t() %>% as.data.frame() %>% as.list()
test <- seqlist[1:2]

hsapiens <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
chimps <- useEnsembl(dataset = "ptroglodytes_gene_ensembl", biomart = "ensembl")
gorilla <- useEnsembl(dataset = "ggorilla_gene_ensembl", biomart = "ensembl")
orangutan <- useEnsembl(dataset = "pabelii_gene_ensembl", biomart = "ensembl")
macaque <- useEnsembl(dataset = "mmulatta_gene_ensembl", biomart = "ensembl")
mouse <- useEnsembl(dataset = "mmusculus_gene_ensembl", biomart = "ensembl")
opossum <- useEnsembl(dataset = "mdomestica_gene_ensembl", biomart = "ensembl")
platypus <- useEnsembl(dataset = "oanatinus_gene_ensembl", biomart = "ensembl")
chicken <- useEnsembl(dataset = "ggallus_gene_ensembl", biomart = "ensembl")

get_sequences_list <- function( ortholist ){
  df <- data.frame(gene_exon_intron = "", ensembl_gene_id = "")
  for(c in 1:length(ortholist)){
    sequ <- getSequence(id = ortholist[[c]], type = "ensembl_gene_id", seqType = "gene_exon_intron", mart = switch(c, 
                                                                                                       "1" = hsapiens,
                                                                                                       "2" = chimps,
                                                                                                       "3" = gorilla,
                                                                                                       "4" = orangutan,
                                                                                                       "5" = macaque,
                                                                                                       "6" = mouse,
                                                                                                       "7" = opossum,
                                                                                                       "8" = platypus,
                                                                                                       "9" = chicken)) %>%
      dplyr::group_by(ensembl_gene_id)
    if(! nrow(sequ) == 0){
      add <- sequ
      df <- dplyr::full_join(df, add)
    }
  }
    df <- df %>% slice(2:n())
    saveRDS(df, file = paste0("Mammal_organs/gene_family_phylogenies/tables/", ortholist[[1]]))
    exportFASTA(df, file = paste0("Mammal_organs/gene_family_phylogenies/fasta_files/", ortholist[[1]], ".fa"))
  
}

mclapply(seqlist, get_sequences_list, mc.cores = 12)

biomartCacheClear()