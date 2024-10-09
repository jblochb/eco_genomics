#vcf filtering 
"coral1"


# Import libraries
{
  library(vcfR)
  library(SNPfiltR)
  library(tidyverse)
  library(qqman)
  library(LEA)
}



setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")
vcf.filt <- hard_filter(vcf, depth = 3)

vcf.filt <- max_depth(vcf.filt, maxdepth = 60)

meta <- read.csv("metadata/meta4vcf.csv", header = T)

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]
names(meta2) <- c("id", "pop")

meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta$id)
options(bitmapType='cairo')
X11.options(type="cairo")
vcf.filt.indMiss <- missing_by_sample(vcf.filt, popmap = meta2, cutoff = 0.45)

vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
#min_mac = minimum minor allele count (this gets rid of NAs)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)
?filter_biallelic
#still need to filter on SNPwise and low freq alleles 
?min_mac
vcf.filt.indSNPmiss <- missing_by_snp(vcf.filt.indMiss, cutoff = 0.5)

DP2 <- extract.gt(vcf.filt.indSNPmiss, element = "DP",
                  as.numeric = T)
str(vcf.filt.indSNPmiss)
heatmap.bp(DP2[1:5000,], rlabels = F, clabels = F)
#DP refers to the overall read depth from all target samples supporting the genotype call 

write.vcf(vcf.filt.indSNPmiss, 
          "~/projects/eco_genomics/population_genomics/outputs/vcf_forty.five_missingness.vcf.gz")

#Diversity Differentiation 
"coral1"



vcf60 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_sixty_missingness.vcf.gz")
vcf45 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_forty.five_missingness.vcf.gz")

meta60 <- meta[meta$id %in% colnames(vcf60@gt[,-1]),]
meta45 <- meta[meta$id %in% colnames(vcf45@gt[,-1]),]
dim(meta45)
{
#calculate the diversity stats using the genetic_diff function in vcfr
#pops = is whatever factor you want to look at diversity across
vcf60.div <- genetic_diff(vcf60,
                        pops = as.factor(meta60$region),
                        method = "nei")
vcf45.div <- genetic_diff(vcf45,
                          pops = as.factor(meta45$region),
                          method = "nei")
?genetic_diff
unique(vcf60.div$CHROM)


chr.main60<- unique(vcf60.div$CHROM)[1:8]
chr.main45<- unique(vcf45.div$CHROM)[1:8]

chrnum60 <- as.data.frame(cbind(chr.main60,seq(1,8,1)))
chrnum45 <- as.data.frame(cbind(chr.main45,seq(1,8,1)))

vcf60.div.MHplot <- left_join(chrnum60, vcf60.div, join_by(chr.main60==CHROM))
vcf60.div.MHplot <- vcf60.div.MHplot %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main60,"_",POS))

vcf60.div.MHplot$V2 = as.numeric(vcf60.div.MHplot$V2)
vcf60.div.MHplot$POS = as.numeric(vcf60.div.MHplot$POS)

vcf45.div.MHplot <- left_join(chrnum45, vcf45.div, join_by(chr.main45==CHROM))
vcf45.div.MHplot <- vcf45.div.MHplot %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(chr.main45,"_",POS))

vcf45.div.MHplot$V2 = as.numeric(vcf45.div.MHplot$V2)
vcf45.div.MHplot$POS = as.numeric(vcf45.div.MHplot$POS)
str(vcf60.div.MHplot)
#manhattan(vcf60.div.MHplot, 
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp = F,
          ylab = "Fst among regions",
          suggestiveline = quantile(vcf60.div.MHplot$Gst, 0.999))

vcf60.div.MHplot %>%
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>%
  filter(value!=0&0.5) %>%
  ggplot(aes(x = value, fill = name))+
  geom_histogram(position="identity", alpha = 0.5, bins = 50)+
  labs(title = "Genome Wide Expected Heterozygosity (Hs)", fill = "Regions",
       x = "Gene Diversity within Regions", y = "Counts of SNPs")
ggsave("Histogram_GenomeDiversity_byRegionofMissingnessSixty.pdf", path = "~/projects/eco_genomics/population_genomics/figures/")

#The following command is used to see avgHs and StdDevHs
#at extremes of this plot where Hst = 0 or 0.5, are likely places where there is only one allele in the population

vcf60.div.MHplot %>%
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>%
  group_by(name) %>%
  filter(value!=0&0.5) %>%
  summarise(avg_Hs=mean(value), StdDev_Hs = sd(value), N_Hs = n())
vcf45.div.MHplot %>%
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>%
  group_by(name) %>%
  filter(value!=0) %>%
  summarise(avg_Hs=mean(value), StdDev_Hs = sd(value), N_Hs = n())
vcf45_Hs = as.data.frame(vcf45.div.MHplot %>%
                           as_tibble() %>% 
                           pivot_longer(c(4:9)) %>%
                           group_by(name) %>%
                           filter(value!=0) %>%
                           summarise(avg_Hs=mean(value), StdDev_Hs = sd(value), N_Hs = n()))
#Next two commands used to see how many SNPs/region have Hs = 0
vcf60_Hs0 = filter(vcf60.div.MHplot, Hs_WEU == 0)
sum(vcf60_Hs0$Hs_WEU)
vcf45_Hs0 = filter(vcf45.div.MHplot, Hs_WEU == 0)
head(vcf45_Hs0)
vcf60.div.MHplot %>%
  select(Hs_CEU,Hs_NE,Hs_NEU,Hs_PNW,Hs_SEU,Hs_WEU)
count(vcf60.div.MHplot, Hs_CEU = 0)
view(vcf45_Hs)
head(vcf60.div.MHplot)
dim(vcf60.div.MHplot)
view(vcf60.div.MHplot)
vcf60_Hs = as.data.frame(vcf60.div.MHplot %>%
                           as_tibble() %>% 
                           pivot_longer(c(4:9)) %>%
                           group_by(name) %>%
                           filter(value!=0&0.5) %>%
                           summarise(avg_Hs=mean(value), StdDev_Hs = sd(value), N_Hs = n()))}
#PCA and Admixture 
"coral1"

vcf.thin60 <- distance_thin(vcf60, min.distance = 500)
vcf.thin45 <- distance_thin(vcf45, min.distance = 500)
?distance_thin
write.vcf(vcf.thin60,"/gpfs1/home/j/b/jblochbe/vcf_final.filtered.thinned.sixty.vcf.gz")
#hide the uncompressed vcf file because it is too big for github (hide outside repo)
write.vcf(vcf.thin60,"~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.sixty.vcf.gz")

write.vcf(vcf.thin45,"/gpfs1/home/j/b/jblochbe/vcf_final.filtered.thinned45.vcf.gz")
#hide the uncompressed vcf file because it is too big for github (hide outside repo)
#write.vcf(vcf.thin60,"~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.thinned.sixty.vcf.gz")

system("gunzip -c ~/vcf_final.filtered.thinned45.vcf.gz > ~/vcf_final.filtered.thinned45.vcf")
setwd("~/projects/eco_genomics/population_genomics/")
geno45 <- vcf2geno(input.file="/gpfs1/home/j/b/jblochbe/vcf_final.filtered.thinned45.vcf",
                   output.file = "outputs/vcf_final.filtered.thinned45.geno")
CentPCA45 <- LEA::pca("outputs/vcf_final.filtered.thinned45.geno", scale = TRUE)
plot(CentPCA45$projections,
     col=as.factor(meta45$region))
legend("bottomright", legend = as.factor(unique(meta45$region)), fill = as.factor(unique(meta45$region)))

CentPCA45<- load.pcaProject(("vcf_final.filtered.thinned45.pcaProject"))
ggplot(as.data.frame(CentPCA45$projections),
       aes(x= V1, y= V2, color = meta45$region, shape= meta45$continent))+
  geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic PCA (45% Individual Missingness Filter)", x = "PC1", y = "PC2", color = "Region", shape = "Color")
# xlim(-10,10)+ ylim(-10,10)
ggsave("figures/CentPCA_PC1vsPC2.45.pdf", width = 6, height = 6, units = "in")

??vcg2geno
system("gunzip -c ~/vcf_final.filtered.thinned.sixty.vcf.gz > ~/vcf_final.filtered.thinned.sixty.vcf")
setwd("~/projects/eco_genomics/population_genomics/")
geno60 <- vcf2geno(input.file="/gpfs1/home/j/b/jblochbe/vcf_final.filtered.thinned.sixty.vcf",
                 output.file = "outputs/vcf_final.filtered.thinned60.geno")
CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned45.geno", scale = TRUE)

plot(CentPCA$projections,
     col=as.factor(meta45$region))
legend("bottomright", legend = as.factor(unique(meta60$region)), fill = as.factor(unique(meta60$region)))
View(CentPCA$projections)
CentPCA<- load.pcaProject(("vcf_final.filtered.thinned45.pcaProject"))
CentPCA$projections$V2 <- -1*CentPCA$projections$V2

CentPCA$projections<- CentPCA$projections %>%
  mutate(V2 = -V2)
?mutate
class(CentPCA$projections)
class(centpca2)
centpca22 <- centpca2 %>% mutate(V2 = -V2)
View(centpca22)
centpca2 <- as.data.frame(CentPCA$projections)

par(mfrow=c(1,3))
#eigenvalues
{PCA75_1<-CentPCA$eigenvalues[1]/sum(CentPCA$eigenvalues)
PCA75_2<-CentPCA$eigenvalues[2]/sum(CentPCA$eigenvalues)

PCA60_1<-CentPCA60$eigenvalues[1]/sum(CentPCA60$eigenvalues)
PCA60_2<-CentPCA60$eigenvalues[2]/sum(CentPCA60$eigenvalues)

PCA45_1<-CentPCA45$eigenvalues[1]/sum(CentPCA45$eigenvalues)
PCA45_2<-CentPCA45$eigenvalues[2]/sum(CentPCA45$eigenvalues)

}
#code block for pca!{
centpca22$eigenvaluse
CentPCA45<- load.pcaProject(("vcf_final.filtered.thinned45.pcaProject"))

p1<-ggplot(as.data.frame(centpca22),
       aes(x= V1, y= V2, color = meta45$region, shape= meta45$continent))+
  geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic PCA (45% Individual Missingness Filter)", x = "PC1 (2.62%)", y = "PC2 (1.25%", color = "Region", shape = "Color")
# xlim(-10,10)+ ylim(-10,
#ggsave("figures/CentPCA_PC1vsPC245.pdf", width = 6, height = 6, units = "in")


#code block for pca60
CentPCA60<- load.pcaProject(("vcf_final.filtered.thinned60.pcaProject"))


p2<- ggplot(as.data.frame(CentPCA60$projections),
       aes(x= V1, y= V2, color = meta60$region, shape= meta60$continent))+
  geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic PCA (60% Individual Missingness Filter)", x = "PC1 (2.42%)", y = "PC2 (1.15%)", color = "Region", shape = "Color")
# xlim(-10,10)+ ylim(-10,10)


CentPCA75<- load.pcaProject(("vcf_final.filtered.thinned.pcaProject"))
p3<- ggplot(as.data.frame(CentPCA75$projections),
       aes(x= V1, y= V2, color = meta2$region, shape= meta2$continent))+
  geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic PCA (75% Individual Missingness Filter)", x = "PC1 (2.30%)", y = "PC2 (1.08%", color = "Region", shape = "Color")
p1 + p2 + p3
pdf("figures/final_pca.pdf")

install.packages("patchwork")
library(patchwork)
?patchwork
"coral1"

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
vcf60 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_sixty_missingness.vcf.gz")
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

meta60 <- meta[meta$id %in% colnames(vcf60@gt[,-1]),]
