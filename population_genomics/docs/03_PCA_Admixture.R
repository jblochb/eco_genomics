library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

?barplot
options(bitmapType = "cairo")
setwd("~/projects/eco_genomics/population_genomics/")
setwd("~/")
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
setwd("~/projects/eco_genomics/population_genomics/")

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
ggsave("figures/CentPCA_PC1vsPC2.75.pdf", width = 6, height = 6, units = "in")
#regraph PC2 vs PC3
ggplot(as.data.frame(CentPCA$projections),
       aes(x= V2, y= V3, color = meta2$region, shape= meta2$continent))+
  geom_point(alpha=0.7)+
  labs(title = "Centaurea genetic PCA", x = "PC2", y = "PC3", color = "Region", shape = "Color")

sum(CentPCA$eigenvalues)
CentPCA$eigenvalues[1]/sum(CentPCA$eigenvalues)
CentPCA$eigenvalues[2]/sum(CentPCA$eigenvalues)
#tells you how much percent a PCA accounts for 


#Admixture analysis aka (Structure) 
#This type of mapping has a genetic model underlying it (HW equilibrium)
#Attempts to groups mixture into # of groups (K) (start with 1, up to maybe 10) which the researcher chooses
#Algorithm assigns each individual to a group, calculate allele freq, calc 2pq, repeat until best result
#Individual might belong to more than 1 group, Q = fractional ancestry

library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)
vcf<- read.vcfR("outputs/vcf_final.filtered.vcf.gz")
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
#must subset the meta final data 

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]

# Now we will run admixture analyses and create admixture plots
#For admixture, we use the LEA R package
#The function inside LEA is called "snmf" (much faster running than Structure)

CentAdmix <- snmf("outputs/vcf_final.filtered.thinned.geno",
                  K = 1:10,
                  entropy = T,
                  repetitions = 3,
                  project = "new",) #if you're adding to this analysis later, you could choose project = "continue"
plot(CentAdmix, color = "blue4", main = "SNMF")
#this plots the cross-entropy score we can use to select models with K values that fit our data well
#investigate where the plot begins to level off but not fully (where the elbow is) in this case: K=3-5
#the relative difference when you add K values is what is significant, you also want to select a less complex model (lower K)
#K = 1:10 means evaluate each K value from 1-10
#entropy = True asks it to test the model 
#The project is saved into :outputs/vcf_final.filtered.thinned.snmfProject 

#To load the project, use:
  #project = load.snmfProject("outputs/vcf_final.filtered.thinned.snmfProject")

#To remove the project, use:
  #remove.snmfProject("outputs/vcf_final.filtered.thinned.snmfProject")
?par
par(mfrow=c(2,1))
plot(CentPCA$eigenvalues[1:10], ylab = "Eigenvalues", xlab = "Number of PCs", color= "blue4", main = "PCA")
plot(CentAdmix, col = "blue4", main = "SNMF")
dev.off()
myK=3
#use to avoid changing in many different places in the future

CE = cross.entropy(CentAdmix, K=myK)
CE
best = which.min(CE)

myKQ = Q(CentAdmix, myK, run = best)
#This Q score is % ancestry in each K group
myKQmeta = cbind(myKQ, meta2)
#cbind can combine two data frames if their rows are identical 
my.colors= c("blue4", "darkolivegreen3", "tomato", "darkorchid", "darkgoldenrod1")

myKQmeta = as_tibble(myKQmeta) %>%
  group_by(continent) %>% 
  arrange(region, pop, .by_group = TRUE)
#first group separately by continent then from within continent by region, within region by group
#tibble is a tidyverse data table
barplot(as.matrix(t(myKQmeta[ , 1:myK])),
        border = NA,
        space = 0,
        col = my.colors[1:myK],
        xlab = "Geographic Regions",
        ylab = "Ancestry Proportions",
        main = paste0("Ancestry Matrix K=",myK))
axis(1,
     at= 1:length(myKQmeta$region), # 1 means mess with x axis instead of y axis
     tick = F,
     labels=myKQmeta$region,
     cex.axis = 0.5, #changes the size of labels
     las = 3) #makes the labels go onto their sides 
pdf("figures/Admixture_K3.pdf")
dev.off()
