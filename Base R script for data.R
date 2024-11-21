## load libraries
library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(patchwork)
library(gridExtra)

options(bitmapType='cairo')
X11.options(type="cairo")

#Load in Data


CentData <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/traits/capitula/PNW_EU_NE_capitulummeasurements.csv")

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/traits/capitula/")
dat <- read.csv("PNW_EU_NE_capitulummeasurements.csv", header=T)

setwd("~/projects/eco_genomics/population_genomics/")
vcf<- read.vcfR("outputs/vcf_final.filtered.vcf.gz")
vcf.thin <- distance_thin(vcf, min.distance = 500)
View(vcf.thin@gt)
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
CentPCA<- load.pcaProject(("vcf_final.filtered.thinned.pcaProject"))
#must subset the meta final data 

###  Data Manipulation and creating PC axes and DF_for_PCA code
{
meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]
dim(meta2)

dat = dat %>%
  mutate(id = paste(Pop, IndID, sep="_"), .before=5) %>%
  group_by(id) %>%
  slice_head(n=1)

pca1 <- prcomp(dat[,c(7:20)], center=T, scale=T)
View(pca1)
pca1.df <- as.data.frame(pca1$x)
pca2 <- prcomp(dat[,c(16,18,19)], center=T, scale=T)
pca2.df <- as.data.frame(pca2$x)


Eigen <- as.data.frame(CentPCA$eigenvalues) # Not sure if this is useful or not
#OR
Projection <- as.data.frame(CentPCA$projections, col = as.factor(meta2$id))
PC1_genomic <- as.data.frame(Projection$V1)
PC1_genomic2 <- as.data.frame(Projection$V2)
PC1_genomic3 <- as.data.frame(Projection$V3)
PC1_trait <- as.data.frame(pca1.df$PC1) #from steve's code
PC1_trait2 <- as.data.frame(pca1.df$PC2)
PC1_trait3 <- as.data.frame(pca1.df$PC3)
PC_trait <- cbind(PC1_trait, PC1_trait2)
PC_trait <- cbind(PC_trait, PC1_trait3)
pc_3var <- as.data.frame(pca2.df$PC1)# Haven't used this to generate second PCA yet
pc_3var2 <- as.data.frame(pca2.df$PC2)
pc_3var3 <- as.data.frame(pca2.df$PC3)
pcs_3var <- cbind(pc_3var, pc_3var2)
pcs_3var <- cbind(pcs_3var, pc_3var3)
View(pca2.df)
 # Now we have pcs_3var, PC_trait, and df_genomic
DF <- as.data.frame(vcf.thin@gt) # Don't think I used this 

DF2 <- as.data.frame(t(DF))
?merge
df_genomic <- cbind(PC1_genomic, PC1_genomic2)
df_genomic <- cbind(df_genomic, PC1_genomic3)
df_genomic <- cbind(meta2, df_genomic)
df_trait <- cbind(PC_trait, pcs_3var)
df_trait <- cbind(dat, df_trait)
df_trait2 <- cbind(dat, pc1_trait2)
DF_for_PCA <- merge(df_genomic, df_trait)
DF_for_PCA2 <- merge(df_genomic, df_trait2)
write.csv(DF_for_PCA, "~/projects/eco_genomics/DF_for_PCA", row.names = F)
# For the first PCA
colnames(DF_for_PCA)[colnames(DF_for_PCA) == 'Projection$V1'] <- 'Genomic_PC1'
colnames(DF_for_PCA)[colnames(DF_for_PCA) == 'Projection$V2'] <- 'Genomic_PC2'
colnames(DF_for_PCA)[colnames(DF_for_PCA) == 'Projection$V3'] <- 'Genomic_PC3'
View(DF_for_PCA)
colnames(DF_for_PCA)[colnames(DF_for_PCA) == 'pca1.df$PC1'] <- 'trait_PC1'
colnames(DF_for_PCA)[colnames(DF_for_PCA) == 'pca1.df$PC2'] <- 'trait_PC2'
colnames(DF_for_PCA)[colnames(DF_for_PCA) == 'pca1.df$PC3'] <- 'trait_PC3'
colnames(DF_for_PCA)[colnames(DF_for_PCA) == 'pca2.df$PC1'] <- 'trait3var_PC1'
colnames(DF_for_PCA)[colnames(DF_for_PCA) == 'pca2.df$PC2'] <- 'trait3var_PC2'
colnames(DF_for_PCA)[colnames(DF_for_PCA) == 'pca2.df$PC3'] <- 'trait3var_PC3'
}
PCA_data <- read.csv("~/projects/eco_genomics/DF_for_PCA")

## No longer needed
{
#For the second PCA (with ratio traits)
colnames(DF_for_PCA2)[colnames(DF_for_PCA2) == 'Projection$V1'] <- 'Genomic_PC1'
View(DF_for_PCA2)
colnames(DF_for_PCA2)[colnames(DF_for_PCA2) == 'pca2.df$PC1'] <- 'trait_PC1'

# For the third version PCA (plotting trait PC2 vs genomic PC1)
colnames(DF_for_PCA2.2)[colnames(DF_for_PCA2.2) == 'Projection$V1'] <- 'Genomic_PC1'
View(DF_for_PCA2.2)
colnames(DF_for_PCA2.2)[colnames(DF_for_PCA2.2) == 'trait_PC1'] <- 'trait_PC2'
}

# Plotting both PCs
PCA1 <- ggplot(PCA_data,
       aes(x = Genomic_PC1, y = trait_PC1, color = region, shape= continent))+
  geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic vs trait PCA", x = "Genomic PC1", y = "Trait PC1", color = "Region", shape = "Color")

#PC1 comparison with 3 variable trait data
PCA2 <- ggplot(PCA_data,
       aes(x = Genomic_PC1, y = trait3var_PC1, color = region, shape= continent))+
  geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic vs trait ratio PCA", x = "Genomic PC1", y = "Trait PC1 (3var)", color = "Region", shape = "Color")

# Plotting PC2 trait vs PC1 Genomics
PCA3 <- ggplot(PCA_data,
       aes(x = Genomic_PC1, y = trait_PC2, color = region, shape= continent))+
  geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic vs trait ratio PCA", x = "Genomic PC1", y = "Trait PC2", color = "Region", shape = "Color")

# Plotting PC2 trait with 3var vs PC1 Genomics
PCA4 <- ggplot(PCA_data,
               aes(x = Genomic_PC1, y = trait3var_PC2, color = region, shape= continent))+
  geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic vs trait ratio PCA", x = "Genomic PC1", y = "Trait PC2 (3var)", color = "Region", shape = "Color")
PCA4

# Plotting PC1 trait vs PC2 Genomics
PCA5 <- ggplot(PCA_data,
               aes(x = Genomic_PC2, y = trait_PC1, color = region, shape= continent))+
  geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic vs trait ratio PCA", x = "Genomic PC2", y = "Trait PC1", color = "Region", shape = "Color")

# Plotting PC1 trait (3var) vs PC2 Genomics
PCA6 <- ggplot(PCA_data,
               aes(x = Genomic_PC2, y = trait3var_PC1, color = region, shape= continent))+
  geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic vs trait ratio PCA", x = "Genomic PC2", y = "Trait PC1 (3var)", color = "Region", shape = "Color")

combined_plot <- grid.arrange(PCA1, PCA2, PCA3, PCA4, PCA5, PCA6, ncol = 2) 

??grid.arrange
#### Probably don't need 
{
PC1_trait2.2 <- as.data.frame(pca2.df$PC2) #from steve's code
df_trait2.2 <- cbind(dat,PC1_trait2.2)

DF_for_PCA2.2 <- merge(df_genomic, df_trait2.2)

trait_PC2 <- 

PCplot <- ggplot(DF_for_PCA2.2,
       aes(x = Genomic_PC1, y = trait_PC2, color = region, shape= continent))+
  geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic vs trait PCA", x = "Genomic PC1", y = "Trait PC2", color = "Region", shape = "Color")
PCplot
}

#################################################
###  For Treemix ###
###################################################

library(vcf2others)
?vcf2treemix
ind_pop <- meta2$region
ind_pop <- factor(ind_pop)
keepers <- factor(unique(ind_pop))

#   #####      Including Missingness
vcf2treemix(
  vcf.thin,
  ind_pop = ind_pop,
  keep_pop = keepers,
  inc_missing = TRUE,
  out_file = "treemix_infile.txt"
)

metatree <- meta[meta$id %in% colnames(vcf.filt3@gt[,-1]),]
ind_pop <- metatree$region
ind_pop <- factor(ind_pop)
keepers <- factor(unique(ind_pop))
vcf2treemix(
  vcf.filt3,
  ind_pop = ind_pop,
  keep_pop = keepers,
  inc_missing = TRUE,
  out_file = "treemix2_infile.txt"
)
treemix <- read.table("treemix_infile.txt")
treemix2 <- read.table("treemix2_infile.txt")
dim("treemix_infile.txt")
summary("treemix_infile.txt")

####### ############ ########### ############ 

##### Investigating filtering for no missingness

##### ########### ######### ############
vcf.filt <- hard_filter(vcf.thin, depth = 3)
vcf.filt <- max_depth(vcf.filt, maxdepth = 40)

vcf.filt.indMiss <- missing_by_sample(vcf.filt, popmap = meta2, cutoff = 0.25)
strwidth(vcf.filt, font = NULL)

vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
#min_mac = minimum minor allele count (this gets rid of NAs)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)

#still need to filter on SNPwise and low freq alleles 

vcf.filt.indSNP <- missing_by_snp(vcf.filt)
vcf.filt.indSNPmiss <- missing_by_snp(vcf.filt.indMiss, cutoff = 0.7)
View(vcf.filt.indSNPmiss)

############ Filter 2 ###########

vcf.filt2 <- missing_by_sample(vcf.filt.indSNPmiss, popmap = meta2, cutoff = 0.20)
vcf.filt2 <- min_mac(vcf.filt2, min.mac = 1)
vcf.filt2 <- missing_by_snp(vcf.filt2, cutoff = 0.8)

vcf.filt3 <- missing_by_sample(vcf.filt2, popmap = meta2, cutoff = 0.15)
vcf.filt3 <- min_mac(vcf.filt3, min.mac = 1)
vcf.filt3 <- missing_by_snp(vcf.filt3, cutoff = 0.9)
dim(vcf.filt3@gt)
#still need to filter on SNPwise and low freq alleles 

individual_missing_proportion <- rowMeans(is.na(vcf.filt.indSNPmiss@gt))

#################

# Imaging the Newick Format Tree

#################

library(ape)
?ape

Phylo_data <- read.tree(text = "(((WEU:0.000703825,((SEU:0.00529387,CEU:0.0016699):0.000857844,PNW:0.00334327):0.00111713):0.00154775,NEU:0.0034095):0.000375701,NE:0.000375701);
")
plot(Phylo_data)
plot.phylo(Phylo_data, type = "unrooted")
Mig_data <- read.tree(text = "((WEU:0.000703824,((SEU:0.0125373,CEU:0.000611748):0.00224212,PNW:0.00399551):0.000464885):0.000773874,(NEU:0.0034095,NE:0.000751401):0.000773874);
0.326489 NA NA NA PNW:0.00399551 SEU:0.0125373")
plot.phylo(Mig_data, type = "unrooted")
setwd("/gpfs1/cl/pbiosw/treemix-1.13/")
# Commands to run on the terminal 
# source("src/plotting_funcs.R")
# plot_tree("treemix_output")