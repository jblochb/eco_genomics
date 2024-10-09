## Code for analyzing RNAseq data using DESeq2

library(DESeq2)
library(ggplot2)
install.packages("DESeq2")
setwd("~/projects/eco_genomics/transcriptomics/")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

BiocManager::install("Biobase", dependencies = T, force = T)

#Import counts matrix 

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = T, row.names = 1)

countsTableRound <- round(countsTable) #Because DESeq2 does not like decimals
tail(countsTable)

#Next is the conditions file 
conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = T, stringsAsFactors = T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~DevTemp + FinalTemp)
