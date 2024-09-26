library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)
options(bitmapType = "cairo")
setwd("~/projects/eco_genomics/population_genomics/")
vcf<- read.vcfR("outputs/vcf_final.filtered.vcf.gz")

# We need to thin the SNPs for LD (linkage disequilibrium) before we run
# PCA and admixture analysis to satisfy the assumptions of independence among loci

vcf.thin <- distance_thin(vcf, min.distance = 500)
#eliminate any SNPs that are closer than 500 BP apart

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
#must subset the meta final data 

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]
dim(meta2)

write.vcf(vcf.thin,"~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz")
#hide the uncompressed vcf file because it is too big for github (hide outside repo)

system("gunzip -c ~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf")

geno <- vcf2geno(input.file = "/gpfs1/home/j/b/jblochbe/vcf_final.filtered.thinned.vcf",
                 output.file = "outputs/vcf_final.filtered.thinned.geno")
CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale = TRUE)
#scale = True makes the mean of each axis = 0 
#LEA::pca - this is a way to tell R which function from which library to use 

plot(CentPCA$projections,
     col=as.factor(meta2$region))
legend("bottomright", legend = as.factor(unique(meta2$region)), fill = as.factor(unique(meta2$region)))

CentPCA<- load.pcaProject(("vcf_final.filtered.thinned.pcaProject"))
show(CentPCA)
plot(CentPCA)
#PC1 is the largest, most dominant eigen value 
#this is a screeplot which shows eigen values in descending order 
ggplot(as.data.frame(CentPCA$projections),
       aes(x= V1, y= V2, color = meta2$region, shape= meta2$continent))+
       geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic PCA", x = "PC1", y = "PC2", color = "Region", shape = "Color")
 # xlim(-10,10)+ ylim(-10,10)
ggsave("figures/CentPCA_PC1vsPC2.pdf", width = 6, height = 6, units = "in")
#regraph PC2 vs PC3
ggplot(as.data.frame(CentPCA$projections),
       aes(x= V2, y= V3, color = meta2$region, shape= meta2$continent))+
  geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic PCA", x = "PC2", y = "PC3", color = "Region", shape = "Color")

sum(CentPCA$eigenvalues)
CentPCA$eigenvalues[1]/sum(CentPCA$eigenvalues)
#tells you how much percent a PCA accounts for 


#Admixture analysis aka (Structure) 
#This type of mapping has a genetic model underlying it (HW equilibrium)
#Attempts to groups mixture into # of groups (K) (start with 1, up to maybe 10) which the researcher chooses
#Algorithm assigns each individual to a group, calculate allele freq, calc 2pq, repeat until best result
#Individual might belong to more than 1 group, Q = fractional ancestry
     