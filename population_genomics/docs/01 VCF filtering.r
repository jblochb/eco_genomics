library(vcfR)

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
getwd()
list.files()

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")

#variants read in the vcf file are SNPs

head(vcf)
#format GT = genotype; AD = read numbers for each allele 

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


