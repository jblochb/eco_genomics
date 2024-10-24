## Code for analyzing RNAseq data using DESeq2

options(bitmapType = "cairo")
X11.options(type="cairo")
{
library(DESeq2)
library(ggplot2)
}#import libraries

install.packages("DESeq2")
setwd("~/projects/eco_genomics/transcriptomics/")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

BiocManager::install("Biobase", dependencies = T, force = T)

#Import counts matrix 

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = T, row.names = 1)

countsTableRound <- round(countsTable) #round to nearest integer because DESeq2 does not like decimals

tail(countsTable)
dim(countsTable)#119439 x 21


#Next is the conditions file 
conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = T, stringsAsFactors = T, row.names = 1)

View(conds)
#########################################################
#Explore Counts Matrix
#########################################################

#Let's see how many reads we have from each sample

colSums(countsTableRound)#outputs the sum of each column which is the total amount of counts/reads for each sample

mean(colSums(countsTableRound)) #18454529 this is a good average for reads/sample for an RNAseq experiment
barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound),
        cex.names = 0.5, las = 2, ylim = c(0,30000000))
abline(h=mean(colSums(countsTableRound)), col = "blue4", lwd=2)

#the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound))#3244.739
median(rowSums(countsTableRound)) #64
#The mean is so much higher than the median
apply(countsTableRound,2,mean)#putting "2" performs this across the columns 
#this is the average amount of reads/transcript for each sample
#gives a sense of variation in sequencing effort across samples

#########################################
#Start analysis in DESeq2

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~DevTemp + FinalTemp)
#consider adding + DevTemp:FinalTemp
dim(dds)

dds <- dds[rowSums(counts(dds) >= 10)>= 15,]
#This is filtering the transcripts, "there must be more than 10 reads for at least 15 samples"
nrow(dds) #35527 = number of transcripts with more than 10 reads for at least 15 samples

#Run the DESeq2 model to test for global differential gene expression
dds <- DESeq(dds)
#DESeq normalizes the data 

#List the results that you've generated
resultsNames(dds)
#[1] "Intercept"             "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28" 
#[4] "FinalTemp_BASE_vs_A28"

#visualize our global gene expression patterns using PCA
#First we need to transform the data for plotting using variance stabilization 

vsd <- vst(dds, blind = F)
#stands for variance stabilized data

pcaData <- plotPCA(vsd, intgroup = cbind("DevTemp", "FinalTemp"), returnData = T)
percentVar <- round(100*attr(pcaData,"percentVar"))#shows variance explained on each axis

final_temp_colors <- c("BASE" = "blue4", "A28" = "hotpink", "A33" = "green4")
shapes_choose <- c("D18" = 16, "D22" = 18)

p <- ggplot(pcaData, aes(PC1, PC2, color = FinalTemp, shape = DevTemp ))+
  geom_point(size = 5)+
  scale_shape_manual(values = shapes_choose)+
  scale_color_manual(values = final_temp_colors)+
  labs(x = paste0("PC1:",percentVar[1], "%"),
       y = paste0("PC2:", percentVar[2], "%"))+
  theme_bw(base_size = 16)
p
