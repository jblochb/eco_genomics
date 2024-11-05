#import libraries
{library(DESeq2)
  library(ggplot2)
  library(tidyverse)
  library(CorLevelPlot)
  library(gridExtra)
  library(eulerr)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(patchwork)
}
setwd("~/projects/eco_genomics/transcriptomics/")
options(bitmapType = "cairo")
X11.options(type="cairo")

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = T, row.names = 1)

countsTableRound <- round(countsTable) #round to nearest integer because DESeq2 does not like decimals

tail(countsTable)
dim(countsTable)#119439 x 21


#Next is the conditions file 
conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = T, stringsAsFactors = T, row.names = 1)

View(conds)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~DevTemp + FinalTemp)
dim(dds)

dds <- dds[rowSums(counts(dds) >= 10)>= 15,]
#This is filtering the transcripts, "there must be more than 10 reads for at least 15 samples"
nrow(dds) #35527 = number of transcripts with more than 10 reads for at least 15 samples

#Run the DESeq2 model to test for global differential gene expression
dds <- DESeq(dds)
#DESeq normalizes the data 

#List the results that you've generated
resultsNames(dds)
#[1] "Intercept"             "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A33" 
#[4] "FinalTemp_BASE_vs_A33"
saveRDS(dds, file = "outputs/dds.rds")
dds <- readRDS("outputs/dds.rds")

res_D22vsD18 <- results(dds, name = "DevTemp_D22_vs_D18", alpha = 0.05)
#order by significance 
res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$padj),]
head(res_D22vsD18)

#look at counts of a specific tip gene that we're interested in to validate that the model is working 
d <- plotCounts(dds,gene = "TRINITY_DN140854_c0_g5_i2", int=(c("DevTemp","FinalTemp")), returnData = T)
d
#here is a way to validate that you're interpreting log2foldchange correctly 
p <- ggplot(d, aes(x= DevTemp, y = count, color = DevTemp, shape = FinalTemp))+
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "gray"))+
  geom_point(position = position_jitter(w=0.2,h=0),size=3)
p

####################################################
#
## Contrasts
#
####################################################

#set up groups within DESeq object
dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
#combining the two factors into one group to aid contrasts we are able to do
design(dds) <- ~group
dds <- DESeq(dds)
dim(dds)
resultsNames(dds)
# [1] "Intercept"               "group_D18A33_vs_D18A28" 
# [3] "group_D18BASE_vs_D18A28" "group_D22A28_vs_D18A28" 
# [5] "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"

#1st CONTRAST! D18 BASE vs D18 A28
{#1. For one DevTemp 18, compare exp at base and A33
res_D18_BASE_D18_A28 <- results(dds, contrast = c("group", "D18BASE", "D18A28"), alpha = 0.05)
res_D18_BASE_D18_A28 <- res_D18_BASE_D18_A28[!is.na(res_D18_BASE_D18_A28$padj),]
res_D18_BASE_D18_A28 <- res_D18_BASE_D18_A28[order(res_D18_BASE_D18_A28$padj),]
head(res_D18_BASE_D18_A28)
summary(res_D18_BASE_D18_A28)

#makek a list of which genes in our comparisons of interest are differentially expressed (List of DEGs)
degs_D18_BASE_D18_A28 <- rownames(res_D18_BASE_D18_A28[res_D18_BASE_D18_A28$padj < 0.05,])

plotMA(res_D18_BASE_D18_A28, ylim=c(1,4))
}
#2nd CONTRAST! D18 BASE vs D18 A33
{#2. For one DevTemp 18, compare exp at base and A33
res_D18_BASE_D18_A33 <- results(dds, contrast = c("group", "D18BASE", "D18A33"), alpha = 0.05)
res_D18_BASE_D18_A33 <- res_D18_BASE_D18_A33[!is.na(res_D18_BASE_D18_A33$padj),]
res_D18_BASE_D18_A33 <- res_D18_BASE_D18_A33[order(res_D18_BASE_D18_A33$padj),]
head(res_D18_BASE_D18_A33)
summary(res_D18_BASE_D18_A33)

#makek a list of which genes in our comparisons of interest are differentially expressed (List of DEGs)
degs_D18_BASE_D18_A33 <- rownames(res_D18_BASE_D18_A33[res_D18_BASE_D18_A33$padj < 0.05,])

plotMA(res_D18_BASE_D18_A33, ylim=c(1,4))}
#3rd CONTRAST! D22 BASE vs D22 A28
{#3. For one DevTemp 22, compare base final temp with A28
  res_D22_BASE_D22_A28 <- results(dds, contrast = c("group", "D22BASE", "D22A28"), alpha = 0.05)
  res_D22_BASE_D22_A28 <- res_D22_BASE_D22_A28[!is.na(res_D22_BASE_D22_A28$padj),]
  res_D22_BASE_D22_A28 <- res_D22_BASE_D22_A28[order(res_D22_BASE_D22_A28$padj),]
  head(res_D22_BASE_D22_A28)
  summary(res_D22_BASE_D22_A28)
  
  #makek a list of which genes in our comparisons of interest are differentially expressed (List of DEGs)
  degs_D22_BASE_D22_A28 <- rownames(res_D22_BASE_D22_A28[res_D22_BASE_D22_A28$padj < 0.05,])
  
plotMA(res_D22_BASE_D22_A28, ylim=c(1,4)) 
}
#4th CONTRAST! D22 BASE vs D22 A33
{#4. For one DevTemp 22, compare base final temp with A33
  res_D22_BASE_D22_A33 <- results(dds, contrast = c("group", "D22BASE", "D22A33"), alpha = 0.05)
  res_D22_BASE_D22_A33 <- res_D22_BASE_D22_A33[!is.na(res_D22_BASE_D22_A33$padj),]
  res_D22_BASE_D22_A33 <- res_D22_BASE_D22_A33[order(res_D22_BASE_D22_A33$padj),]
  head(res_D22_BASE_D22_A33)
  summary(res_D22_BASE_D22_A33)
  
  #makek a list of which genes in our comparisons of interest are differentially expressed (List of DEGs)
  degs_D22_BASE_D22_A33 <- rownames(res_D22_BASE_D22_A33[res_D22_BASE_D22_A33$padj < 0.05,])
  
  plotMA(res_D22_BASE_D22_A33, ylim=c(1,4)) 
}

#to see how many differential expressed genes for each contrast
length(degs_D18_BASE_D18_A28) #41
length(degs_D18_BASE_D18_A33) #332

length(degs_D22_BASE_D22_A28) #289
length(degs_D22_BASE_D22_A33) #1564

length(intersect(degs_D18_BASE_D18_A28,degs_D18_BASE_D18_A33)) #34
41 - 34 #7 DEGs between base and A28 at 18 DevTemp
332 - 34 #298 DEGs between base and A33 at 18 DevTemp

Euler18 <- euler(c("A28"=7,"A33"=298,"A28&A33"=34))
Euler18
Euler1 <- plot(Euler18, lty=1:3, quantities =T)

length(intersect(degs_D22_BASE_D22_A28,degs_D22_BASE_D22_A33)) #144
289 - 144 #145 DEGs between base and A28 at 22 DevTemp
1564 - 144 #1420 DEGs between base and A33 at 22 DevTemp
D22diffs <- as_data_frame(intersect(degs_D22_BASE_D22_A28,degs_D22_BASE_D22_A33))
Euler22 <- euler(c("A28"=145,"A33"=1420,"A28&A33"=144))
Euler22
Euler2 <- plot(Euler22, lty=1:3, quantities = T)
?eulerr
?grid.arrange
combined_plot <- grid.arrange(Euler1, Euler2, ncol = 2)
#################################
#
# Make a scatter plot of responses for each devtemp
#
#################################

#Scatter plot for D18
{
# contrast D18_A28vsBASE
res_D18_BASEvsA28 <- as.data.frame(results(dds, contrast= c("group","D18BASE","D18A28"), alpha = 0.05))

# contrast D18_A33vsBASE
res_D18_BASEvsA33 <- as.data.frame(results(dds, contrast= c("group","D18BASE","D18A33"), alpha = 0.05))

#Merge dataframes

res_dfBASE <- merge(res_D18_BASEvsA28, res_D18_BASEvsA33, by = "row.names", suffixes = c(".A28",".A33"))
rownames(res_dfBASE) <- res_dfBASE$Row.names
res_dfBASE <- res_dfBASE[,-1]

res_dfBASE <- res_dfBASE %>%
  mutate(fill= case_when(
    padj.A28<0.05 & stat.A28 <0 ~ "tomato", # genes that are significant and upregulated in experimental env
    padj.A28<0.05 & stat.A28 >0 ~ "darkred",
    padj.A33<0.05 & stat.A33 <0 ~ "cadetblue",
    padj.A33<0.05 & stat.A33 >0 ~"cadetblue1"
  ))

# plot above

#Count the number of points per fill color 

color_counts <- res_dfBASE %>% 
  group_by(fill) %>%
  dplyr::summarize(count = dplyr::n())
label_positions <- data.frame(
  fill = c("cadetblue", "darkred", "cadetblue1","tomato"),
  x_pos =c(1,5,0,-7.5),
  y_pos =c(-5,0,9,3)
)

label_data <- merge(color_counts, label_positions, by = "fill")
plotBASE <- ggplot(res_dfBASE, aes(x = log2FoldChange.A28, y = log2FoldChange.A33, color = fill))+
  geom_point(alpha = 0.8)+
  scale_color_identity()+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey")+
  xlim(-10,10)+ ylim(-10,10)+
  labs(x = "LFC (BASE vs 28°C at 18°C)", y = "LFC (BASE vs 33°C at 18°C)")+
  theme_minimal()+
  geom_text(data = label_data, aes(x = x_pos, y = y_pos, label = count, color = fill), size = 5)
plotBASE
}

#Scatter plot for D22
{
# contrast D18_A28vsBASE
res_D22_BASEvsA28 <- as.data.frame(results(dds, contrast= c("group","D22BASE","D22A28"), alpha = 0.05))

# contrast D18_A33vsBASE
res_D22_BASEvsA33 <- as.data.frame(results(dds, contrast= c("group","D22BASE","D22A33"), alpha = 0.05))

#Merge dataframes

res_dfBASE22 <- merge(res_D22_BASEvsA28, res_D22_BASEvsA33, by = "row.names", suffixes = c(".A28",".A33"))
rownames(res_dfBASE22) <- res_dfBASE22$Row.names
res_dfBASE22 <- res_dfBASE22[,-1]

res_dfBASE22 <- res_dfBASE22 %>%
  mutate(fill= case_when(
    padj.A28<0.05 & stat.A28 <0 ~ "tomato", # genes that are significant and upregulated in experimental env
    padj.A28<0.05 & stat.A28 >0 ~ "darkred",
    padj.A33<0.05 & stat.A33 <0 ~ "cadetblue",
    padj.A33<0.05 & stat.A33 >0 ~"cadetblue1"
  ))

# plot above

#Count the number of points per fill color 

color_counts22 <- res_dfBASE22 %>% 
  group_by(fill) %>%
  dplyr::summarize(count = dplyr::n())
label_positions <- data.frame(
  fill = c("cadetblue", "darkred", "cadetblue1","tomato"),
  x_pos =c(2,5,0,-5),
  y_pos =c(-6,0,9,1)
)

label_data22 <- merge(color_counts22, label_positions, by = "fill")
plotBASE22 <- ggplot(res_dfBASE22, aes(x = log2FoldChange.A28, y = log2FoldChange.A33, color = fill))+
  geom_point(alpha = 0.8)+
  scale_color_identity()+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "grey")+
  xlim(-10,10)+ ylim(-10,10)+
  labs(x = "LFC (BASE vs 28°C at 22°C)", y = "LFC (BASE vs 33°C at 22°C)")+
  theme_minimal()+
  geom_text(data = label_data22, aes(x = x_pos, y = y_pos, label = count, color = fill), size = 5)
plotBASE22
}
combined_plot <- grid.arrange(Euler1, Euler2, ncol = 2) 
patchwork<- plotBASE + plotBASE22

patchwork + plot_annotation(
  tag_levels = "A"
)
pdf("figures/final_LFCplot.pdf")
ggsave("~/projects/eco_genomics/transcriptomics/figures/combined_scatter_plotBASE18and22.png", combined_plot, width = 12, height = 6)
help("patchwork")

?ggplot

#Investigating Base18 vs Base22

res_D18_BASE_D22_BASE <- results(dds, contrast = c("group", "D18BASE", "D22BASE"), alpha = 0.05)
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[!is.na(res_D18_BASE_D22_BASE$padj),]
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[order(res_D18_BASE_D22_BASE$padj),]
head(res_D18_BASE_D22_BASE)
summary(res_D18_BASE_D22_BASE)
degs_D18_BASE_D22_BASE <- rownames(res_D18_BASE_D22_BASE[res_D18_BASE_D22_BASE$padj < 0.05,])
length(degs_D18_BASE_D22_BASE)

citation(package = "gridExtra", lib.loc = NULL)
