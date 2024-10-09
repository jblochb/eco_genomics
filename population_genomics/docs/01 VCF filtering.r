library(vcfR)

#set working directory to where the data file lives 
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
getwd()

#you can put any folder into the paratheses of this command
list.files()

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")
vcf
#variants read in the vcf file are SNPs

head(vcf)
#format GT = genotype; PL = a phred score; DP = the depth; AD = read numbers for each allele 
#when there is no data it's coded as ./.

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format = "fasta")

gff<- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep = "\t", quote = "")

chr1<- create.chromR(name="Chromosome 1", vcf = vcf, seq = dna, ann=gff)

plot(chr1)

chromoqc(chr1, xlim = c(1e1, 1.1e8))
chromoqc(chr1, xlim = c(1e5, 1.1e6))
#red bars show where the genes are in the genome 

#use below to save to project direct
pdf(file="~/projects/eco_genomics/population_genomics/figures/ChromoPlot.pdf")
chromoqc(chr1, xlim = c(1e1, 1.1e8))
dev.off()

DP <- extract.gt(vcf, "DP", as.numeric = T)
dim(DP)
#rows is the number of SNPs and columns is the number of individuals 
DP[1:5,1:10]
quantile(DP)

DP[DP==0] <- NA
#for where locations inside DP are 0, rename them as NA
quantile(DP, na.rm = T)

# visualize the matrix of DP and missingness in our VCF file:
?heatmap.bp
heatmap.bp(DP)
heatmap.bp(DP[1:1000,], rlabels = F, clabels = F)

library(SNPfiltR)

?SNPfiltR
#exploring filtering by depth
hard_filter(vcf)
vcf.filt<- vcf
vcf.filt <- hard_filter(vcf, depth = 3)

vcf.filt <- max_depth(vcf.filt, maxdepth = 60)
#displays that the mean depth for SNPs is about 30. Adding the 60 filters out genotypes with greater than 60 reads

meta <- read.csv("metadata/meta4vcf.csv", header = T)
View(meta)
#names in the meta file must be exactly the same as the vcf file 
meta2 <- meta[,c(1,4)]
#above command is subsetting the data. There is nothing before the comma because we want all the rows
head(meta2)
#rename region to pop to use the SNP function 
names(meta2) <- c("id", "pop")

meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta$id)
options(bitmapType='cairo')
?X11.options
X11.options(type="cairo")
vcf.filt.indMiss <- missing_by_sample(vcf.filt, popmap = meta2, cutoff = 0.60)
strwidth(vcf.filt, font = NULL)
?missing_by_sample
?strwidth
dim(vcf.filt)
vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
#min_mac = minimum minor allele count (this gets rid of NAs)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)

#still need to filter on SNPwise and low freq alleles 

vcf.filt.indSNP <- missing_by_snp(vcf.filt)
vcf.filt.indSNPmiss <- missing_by_snp(vcf.filt.indMiss, cutoff = 0.5)

DP2 <- extract.gt(vcf.filt.indSNPmiss, element = "DP",
                  as.numeric = T)
heatmap
heatmap.bp(DP2[1:5000,], rlabels = F, clabels = F)

write.vcf(vcf.filt.indSNPmiss, 
           "~/projects/eco_genomics/population_genomics/outputs/vcf_final.filtered.vcf.gz")

#  ~/ is a shortcut that means home directory 

#filter_biallelic removes multi-allelic sites
vcf.filt.indMiss <- filter_biallelic(vcf.filt)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)
vcf.filt.indSNPmiss <- missing_by_snp(vcf.filt.indMiss, cutoff = 0.5)
