CentData <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/traits/capitula/PNW_EU_NE_capitulummeasurements.csv")


setwd("/gpfs1/cl/pbio3990/PopulationGenomics/traits/")

View(CentData)

setwd("~/projects/eco_genomics/population_genomics/")
vcf<- read.vcfR("outputs/vcf_final.filtered.vcf.gz")
vcf.thin <- distance_thin(vcf, min.distance = 500)
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
#must subset the meta final data 

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]
dim(meta2)

Eigen <- as.data.frame(CentPCA$eigenvalues)
#OR
PC1_genomic <- as.data.frame(Projection$V1)
PC1_trait <- as.data.frame(pca1.df$PC1) #from steve's code
pca2.df <- as.data.frame(pca2$x)
pc1_trait2 <- as.data.frame(pca2.df$PC1)# Haven't used this to generate second PCA yet


??pcaProject
Projection <- as.data.frame(CentPCA$projections, col = as.factor(meta2$id))
DF <- as.data.frame(vcf.thin@gt)

DF2 <- as.data.frame(t(DF))
?merge
df_genomic <- cbind(meta2, PC1_genomic)
df_trait <- cbind(dat,PC1_trait)
DF_for_PCA <- merge(df_genomic, df_trait)

colnames(DF_for_PCA)[colnames(DF_for_PCA) == 'Projection$V1'] <- 'Genomic_PC1'
View(DF_for_PCA)
colnames(DF_for_PCA)[colnames(DF_for_PCA) == 'pca1.df$PC1'] <- 'trait_PC1'

ggplot(DF_for_PCA,
       aes(x = Genomic_PC1, y = trait_PC1, color = region, shape= continent))+
  geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic vs trait PCA", x = "Genomic PC1", y = "Trait PC1", color = "Region", shape = "Color")
?cbind

