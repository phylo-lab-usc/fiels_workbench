#!/usr/bin/env Rscript

# Custom code - Script A : 
# Gene expression dynamics during rapid organismal diversification of East African cichlid fishes.            
# Normalisation and count matrix 

# Import DESEQ2:
library("DESeq2")

# Path to count files :
Path_toCount = '/scicore/home/salzburg/eltaher/cichlidX/data/ReadCount'

## Define the condition and the individus:
count_name = list.files(path = Path_toCount,full.names=F)
sample_name = unlist(lapply(count_name,function(x) strsplit(x,'.count')[[1]][1]))


# SPECIMEN INFORMATION: 
# 1. Species id
# 2. Tissu type
# 3. Sex

sp = unlist(lapply(sample_name,function(x) strsplit(x,'_')[[1]][1]))
tissu = unlist(lapply(sample_name,function(x) strsplit(x,'_')[[1]][4]))
sex = unlist(lapply(sample_name,function(x) strsplit(x,'_')[[1]][5]))

# Information Table: 
sampleTable = data.frame(
  sampleName = sample_name,
  fileName = count_name,
  condition_sp = sp,
  condition_tissu = tissu,
  condition_sex = sex)

## Function that generate the count_matrix:
ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory = Path_toCount ,design= ~ condition_sp + condition_tissu + condition_sex )
ddsHTSeq_new = ddsHTSeq

# Remove lowely expressed genes (less than 5 reads in 3 samples, DESEQ2 recommendation):
ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)>5) >= 3, ]

# VST normalisation: 
rld <- vst(ddsHTSeq , blind=FALSE)
norm_count = data.frame('contigName' = rownames(assay(rld)))
norm_count = cbind(norm_count,assay(rld))
dds_Wald = DESeq(ddsHTSeq, test="Wald")
# Save working space (vst normalization takes some time)
save.image('/scicore/home/salzburg/eltaher/cichlidX/result/DESEQ2/CichlidX_GeneExpressionEvolution_script_A.Rdata')
