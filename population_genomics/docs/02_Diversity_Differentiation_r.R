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
#allelic divergence between populations (also the amount of non random mating that's occuring')
manhattan(vcf.div.MHplot, 
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("blue4","orange3"),
          logp = F,
          ylab = "Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.999))
