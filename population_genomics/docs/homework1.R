#vcf filtering 
"coral1"

library(vcfR)
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")
library(SNPfiltR)
vcf.filt <- hard_filter(vcf, depth = 3)

vcf.filt <- max_depth(vcf.filt, maxdepth = 60)

meta <- read.csv("metadata/meta4vcf.csv", header = T)
meta2 <- meta[,c(1,4)]
names(meta2) <- c("id", "pop")

meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta$id)
options(bitmapType='cairo')
X11.options(type="cairo")
vcf.filt.indMiss <- missing_by_sample(vcf.filt, popmap = meta2, cutoff = 0.45)

vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
#min_mac = minimum minor allele count (this gets rid of NAs)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)

#still need to filter on SNPwise and low freq alleles 

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

library(vcfR)
library(tidyverse)
library(qqman)

vcf60 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_sixty_missingness.vcf.gz")
vcf45 <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_forty.five_missingness.vcf.gz")

meta60 <- meta[meta$id %in% colnames(vcf60@gt[,-1]),]
meta45 <- meta[meta$id %in% colnames(vcf45@gt[,-1]),]
dim(meta45)

#calculate the diversity stats using the genetic_diff function in vcfr
#pops = is whatever factor you want to look at diversity across
vcf60.div <- genetic_diff(vcf60,
                        pops = as.factor(meta60$region),
                        method = "nei")
vcf45.div <- genetic_diff(vcf45,
                          pops = as.factor(meta45$region),
                          method = "nei")

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
                           summarise(avg_Hs=mean(value), StdDev_Hs = sd(value), N_Hs = n()))
#PCA and Admixture 
"coral1"

library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

vcf.thin60 <- distance_thin(vcf60, min.distance = 500)
vcf.thin45 <- distance_thin(vcf45, min.distance = 500)
