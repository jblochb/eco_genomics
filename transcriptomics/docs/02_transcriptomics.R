#load in day 1 transcriptomic script data and run DEseq object
library(pheatmap)

#List the results that you've generated 
resultsNames(dds)
#[1] "Intercept"             "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28" 
#[4] "FinalTemp_BASE_vs_A28"
saveRDS(dds, file = "outputs/dds.rds")
#Pull out the results for developmental temp 22 vs 18 degrees
res_D22vsD18 <- results(dds, name = "DevTemp_D22_vs_D18", alpha = 0.05)
#order by significance 
res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$padj),]
head(res_D22vsD18)
# the first transcript in this object has the highest difference in expression level
#within the result we chose 
#log2foldchange is a measure of difference in gene expression 
#devtemp 22 has about 2x less expression than devtemp 18

#look at counts of a specific tip gene that we're interested in to validate that the model is working 
d <- plotCounts(dds,gene = "TRINITY_DN140854_c0_g5_i2", int=(c("DevTemp","FinalTemp")), returnData = T)
d

#here is a way to validate that you're interpreting log2foldchange correctly 
p <- ggplot(d, aes(x= DevTemp, y = count, color = DevTemp, shape = FinalTemp))+
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "gray"))+
  geom_point(position = position_jitter(w=0.2,h=0),size=3)
p

#MA Plot 
plotMA(res_D22vsD18, ylim=c(-4,4))
#shows that some genes that have high expression (house keeping genes) usually have
#less significant log fold change 

#volcano plot

#convert our deseq results object into a data frame to plot
res_df <- as.data.frame(res_D22vsD18)
#then we'll add a column to the data frame to say if gene is significantly differentially expressed or not
res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >1, "Significant", "Not Significant")
View(res_df)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significant))+
  geom_point(alpha = 0.8)+
  scale_color_manual(values = c("slateblue", "tomato"))+
  labs(x = "Log2 Fold Change", y = "Log10 Adjusted P-value", title = "Volcano Plot")+
  theme_minimal()+
  theme(legend.position = "top")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "orange")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", colour = "orange")

#########           Heat  Map

vsd <- vst(dds, blind=F)

topgenes <- head(rownames(res_D22vsD18),20)
mat <- assay(vsd)[topgenes,]
df <- as.data.frame(colData(dds)[,c("DevTemp","FinalTemp")])
pheatmap(mat, annotation_col = df, show_rownames = F, cluster_cols=T,cluster_rows=T)
