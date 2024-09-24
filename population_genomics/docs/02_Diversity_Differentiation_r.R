#Estimating diversity and genetic differentiation in the filtered Centauria data

library(vcfR)
library(tidyverse)
library(qqman)

# helps solve plotting issues 
X11.options(type="cairo")

#read in our VCF file from our repo outputs directory 

vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

#read in our metadata - - Info on population of origin, what regions pops come from, wat continent

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
head(meta)
vcf
#vcf file has 595 samples 
dim(meta)
#this command gives the dimensions produces 629 individuals; we need to filter this to make it compatible with vcf

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]
# %in% is an operator that says find all of the values in common and exports the identical values into a new dataframe

vcf@gt[1:5,1:3]
vcf@gt[1:3,1:5]
#@ svmbol can be used to navigate to the headings of dataframes within a larger file
colnames(vcf@gt[,-1])
colnames(vcf@gt[,-1])[1]
colnames(vcf@gt[,])
dim(meta2)

#calculate the diversity stats using the genetic_diff function in vcfr
#pops = is whatever factor you want to look at diversity across
vcf.div <- genetic_diff(vcf,
                        pops = as.factor(meta2$region),
                        method = "nei")

str(vcf.div)
#allows you to see the structure of the file

unique(vcf.div$CHROM)
#lists the different values/words
#in the reference genome, chroms that start with cm are real chromosomes
#chromosomes that start with JARY are scaffolds, unsure where these strands go
chr.main<- unique(vcf.div$CHROM)[1:8]
#use the colon to mean 1-8

#has to tell qqman each chromosome number
chrnum <- as.data.frame(cbind(chr.main,seq(1,8,1)))
#seq is used to rename words to numbers 
#try to rename v-2 column?
#tidy operation that joins two dataframe
vcf.div.MHplot <- left_join(chrnum, vcf.div, join_by(chr.main==CHROM))
head(vcf.div.MHplot)
vcf.div.MHplot <- vcf.div.MHplot %>%
        filter(Gst>0) %>%
        mutate(SNP=paste0(chr.main,"_",POS))
#some of the variables are not numerical (R is reading them as characters)
vcf.div.MHplot$V2 = as.numeric(vcf.div.MHplot$V2)
vcf.div.MHplot$POS = as.numeric(vcf.div.MHplot$POS)
str(vcf.div.MHplot)
manhattan(vcf.div.MHplot, 
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp = F,
          ylab = "Fst among regions")
#Each dot is a SNP, colors alternate between chromosomes, v axis = Fst which means 
#allelic divergence between populations (also the amount of non random mating that's occurring')
manhattan(vcf.div.MHplot, 
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp = F,
          ylab = "Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.999))
#some SNPs stand out in terms of genetic differentiation, they are likely under selection

write.csv(vcf.div.MHplot, "~/projects/eco_genomics/population_genomics/outputs/Genetic_Diff_By_Region.csv",
          quote = F,
          row.names = F)
#where you want the file to be in quotes and add the name you want
#no row names and no quotes to prevent future errors

names(vcf.div.MHplot)
#Hs values are stored in columns 4-9

#We want to create one column that contains all the Hs values from each column
vcf.div.MHplot %>%
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x = value, fill = name))+
  geom_histogram(position="identity", alpha = 0.5, bins = 50)+
  labs(title = "Genome Wide Expected Heterozygosity (Hs)", fill = "Regions",
       x = "Gene Diversity within Regions", y = "Counts of SNPs")
ggsave("Histogram_GenomeDiversity_byRegion.pdf", path = "~/projects/eco_genomics/population_genomics/figures/")

#at extremes of this plot where Hst = 0 or 0.5, are likely places where there is only one allele in the population

vcf.div.MHplot %>%
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>%
  group_by(name) %>%
  filter(value!=0&0.5) %>%
  summarise(avg_Hs=mean(value), StdDev_Hs = sd(value), N_Hs = n())
#summarise can be used for many different statistics
#with filtering out the values of 0, we can see the average Hs when we exclude the two extremes 
#In R you can use != for does not equal 
X11.options(type="cairo")
options(bitmapType = "cairo")
#Imaging the scaffolds 

library(vcfR)
library(tidyverse)
library(qqman)

vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

unique(vcf.div$CHROM)
scaffolds<- unique(vcf.div$CHROM)[9:17]

#has to tell qqman each chromosome number
scaffnum <- as.data.frame(cbind(scaffolds,seq(1,9,1)))
view(scaffnum)
vcf.div.MHplot <- left_join(scaffnum, vcf.div, join_by(scaffolds==CHROM))

vcf.div.MHplot <- vcf.div.MHplot %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(scaffolds,"_",POS))

vcf.div.MHplot$V2 = as.numeric(vcf.div.MHplot$V2)
vcf.div.MHplot$POS = as.numeric(vcf.div.MHplot$POS)
str(vcf.div.MHplot)
view(vcf.div.MHplot)
manhattan(vcf.div.MHplot, 
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp = F,
          ylab = "Fst among scaffold regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.999))
head(vcf.div.MHplot)
?head
head(vcf.div.MHplot, n = 15)
print(vcf.div.MHplot$V2)
manhattan(vcf.div.MHplot, 
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp = F,
          ylab = "Fst among scaffold regions")
print(vcf.div.MHplot$POS)
?manhattan
#maybe the POS are very different for each scaffold because they do not map well onto the manhattan plot 

#image the scaffolds alongside the chromosomes?

library(vcfR)
library(tidyverse)
library(qqman)

vcf <- read.vcfR("~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

unique(vcf.div$CHROM)
chrm.name<- unique(vcf.div$CHROM)[1:17]

#has to tell qqman each chromosome number
chrmnum <- as.data.frame(cbind(chrm.name,seq(1,17,1)))
view(chrmnum)
vcf.div.MHplot <- left_join(chrmnum, vcf.div, join_by(chrm.name==CHROM))

vcf.div.MHplot <- vcf.div.MHplot %>%
  filter(Gst>0) %>%
  mutate(SNP=paste0(chrm.name,"_",POS))

vcf.div.MHplot$V2 = as.numeric(vcf.div.MHplot$V2)
vcf.div.MHplot$POS = as.numeric(vcf.div.MHplot$POS)
str(vcf.div.MHplot)

manhattan(vcf.div.MHplot, 
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp = F,
          ylab = "Fst among scaffold regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.999))
#still difficult to make out