## Script for analyzing and visualizing gene correlation networks 

#import libraries
{library(DESeq2)
library(ggplot2)
library(WGCNA); options(stringsAsFactors = F)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(Rmisc)
  }
setwd("~/projects/eco_genomics/transcriptomics/")
options(bitmapType = "cairo")
X11.options(type="cairo")

#Countstable and conds copied from 01_DESeq file
# STEP 1 importing our data
{countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = T, row.names = 1)

countsTableRound <- round(countsTable) #round to nearest integer because DESeq2 does not like decimals

tail(countsTable)
dim(countsTable)#119439 x 21


#Next is the conditions file 
conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = T, stringsAsFactors = T, row.names = 1)}

traitData = read.table("/gpfs1/cl/pbio3990/Trait_Data.txt", header = T, row.names = 1)

#Filter the matrix to just BASE data (because those are the data for which we have traits measured )

filtered_count_matrix_BASEonly <- countsTable[,conds$FinalTemp == "BASE"]
filtered_sample_metadata_BASEonly <- conds[conds$FinalTemp== "BASE",]
rounded_filtered_count_matrix <- round(filtered_count_matrix_BASEonly)

# STEP 2 Detecting outliers 
# detecting outlier genes 
gsg <- goodSamplesGenes(t(rounded_filtered_count_matrix))
summary(gsg)
table(gsg$goodGenes) # some of the initial transcripts that we imported in may have over-dispersion or none or low expression 
#this flags the genes that are not suitable for this analysis 
# FALSE  TRUE 
# 37235 82203 
table(gsg$goodSamples) # All are good! 

#Filter out the bad genes 
data_WGCNA <- rounded_filtered_count_matrix[gsg$goodGenes == T,]
dim(data_WGCNA)

# Use clustering with a tree dendrogram to identify outlier samples 
htree <- hclust(dist(t(data_WGCNA)), method = "average")
plot(htree)

# PCA- another outlier detection method 

pca <- prcomp(t(data_WGCNA))
pca_data <- pca$x
#transform it into a data frame 
pca_data <- as.data.frame(pca_data)

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
ggplot(pca_data, aes(PC1,PC2))+
  geom_point()+
  geom_text(label = rownames(pca_data))+
  labs(x = paste0("PC1: ", pca.var.percent[1], " %"), y = paste0("PC2: ", pca.var.percent[2]," %"))

# STEP 3 - Normalization

colData <- row.names(filtered_sample_metadata_BASEonly)

dds_WGCNA <- DESeqDataSetFromMatrix(countData = data_WGCNA,
                                    colData = filtered_sample_metadata_BASEonly,
                                    design = ~1) #there are no specified groups 
dds_WGCNA_75 <- dds_WGCNA[rowSums(counts(dds_WGCNA) >= 15) >= 6,]
#the sum of the counts per row (which I think each row is a transcript) must be more than or equal to 15 for at least 6 of the samples 
nrow(dds_WGCNA_75) #filtered down to 29,559 transcripts

dds_norm <- vst(dds_WGCNA_75) # We are performing the variance stabilization function on our data 

# get and save normalized counts to use below 
norm.counts <- assay(dds_norm) %>% 
  t()
# STEP 4- Network Connection

# Choose a set of soft thresholding powers 
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

#Call the network topology analysis function (takes a couple minutes to run)
sft <- pickSoftThreshold(norm.counts, powerVector = power,
                         networkType = "signed",
                         verbose = 5)
# "signed" means only paying attention to genes that are positively correlated with each other and not negatively related (both up or both down)

sft.data <- sft$fitIndices

# plot to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept = 0.8, color = "red")+
  labs(x = "Power", y = "Scale free topology model fit, signed R^2")+
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept = 0.8, color = "red")+
  labs(x = "Power", y = "Mean Connectivity")+
  theme_classic()

grid.arrange(a1, a2, nrow = 2)
# we will choose a power value of 26

soft_power <- 26
temp_cor <- cor
cor <- WGCNA::cor # this sets the temp_cor function to use WGCNA's correlation function 

norm.counts[] <- sapply(norm.counts, as.numeric) #creates a numerical data table

# The command below creates the network and identifies modules based on the parameters that we chose 
bwnet26 <- blockwiseModules(norm.counts, maxBlockSize = 30000,
                            TOMType = "signed",
                            power = soft_power,
                            mergeCutHeight = 0.25,
                            numericLabels = F,
                            randomSeed = 1234,
                            verbose = 3) #whether you allow correlations to be only positive or both signed correlations 
# in this case TOMType will only focus on positive correlations because these are more biologically relevant

saveRDS(bwnet26, file = "outputs/bwnet26.rds")
# To load the bwnet file in later use:
bwnet26 <- readRDS("outputs/bwnet26.rds") # 

cor <- temp_cor # This resets the cor function to base R's cor function instead of using WGCNA's cor function

#STEP FIVE: Explore module eigengenes 

module_eigengenes <- bwnet26$MEs
head(module_eigengenes)
dim(modul_eigengenes)

#Get the number of genes for each module 

table(bwnet26$colors) #because it names each of the modules with a random color

#Plot the dendrogram and the module colors 
plotDendroAndColors(bwnet26$dendrograms[[1]], cbind(bwnet26$unmergedColors, bwnest26$colors),
                    c("unmerged", "merged"),
                    dendroLabels = F,
                    addGuide = T,
                    hang = 0.03,
                    guideHang = 0.05)

# Step SIX: Correlation of modules with traits 

# Define the numbers of genes and samples 
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

# This will test for a correlation between module eigengenes and trait data
module.trait.corr <- cor(module_eigengenes, traitData, use = "p")

# Calculate p values for each correlation 
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# Visualize module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traitData, by = "row.names")
head(heatmap.data) # last three columns are trait data 

# Address error of row.names not being numeric
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = "Row.names")
names(heatmap.data)

# Make pretty heatmap of correlations 
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[42:44], # these values may need to change based on 
             y = names(heatmap.data)[1:41],
             col = c("blue2", "skyblue", "white", "pink", "red"))# number of eigengenes
