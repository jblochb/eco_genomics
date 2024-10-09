BiocManager:: install("DESeq2", dependencies = T, force = T)

setwd("~/projects/eco_genomics/transcriptomics/")

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt")

countsTableRound <- round(countsTable)

