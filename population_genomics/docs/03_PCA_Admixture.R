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

meta2 <- meta[meta$id %in% colnames(vcf.thin@gt[,-1]),]
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

     