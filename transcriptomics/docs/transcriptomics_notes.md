# Coding and data notes for the Transcriptomics module

## Author: John Blochberger

### 10/8/2024 - Transcriptomics intro and DESeq2

We began with a lecture on transcriptomics and background info about our dataset. In our data set there are two factors with their respective levels.

-   Development Temp

    -   18 degrees Celsius and 22 degrees Celsius

-   Final Temp

    -   Baseline temp, 28 degrees Celsius, and 33 degrees Celsius

We troubleshooted downloading the DESeq2 package.

### 10/10/2024 - Transcriptomics day 2

Today we started with the counts matrix/table. My notebook has a good diagram that can help me understand the layout of this file. The columns of this file are divided by each sample and its replicates. For this study, each column is a combination of devTemp and FinalTemp. Each row is a different transcript and the values are the amount of reads that mapped to each transcript. 
We found the total amount of reads per column (sample) and averaged them to find an average of 18,454,529 reads per sample.This was deemed a good amount of reads for an RNAseq experiment. We made a bar plot so we could visually compare how each sample's data compared to the mean and to each other. If we do the mean of the row sums then we get the average amount of reads for each transcript. In this case it was about 3,200, and the median was only about 64. This shows that we have massive variance and non-normally distributed data which is typical for RNAseq. 

We then used DESeqDataSetFromMatrix to combine the conditions data with the rounded counts table to create a dataset for DESeq to work on. We applied a filter to this to remove very lowly expressed transcripts. The result was 35527 (number of transcripts with more than 10 reads for at least 15 samples). We used vst() create "variance stabilized data" and imaged this on a PCA. We could clearly see that the different final temperature samples clustered in PC space. 


### 10/15/2024 - Transcriptomics day 3 ()

We created plots that imaged the log2foldchange so that we could be sure of our reference point and that we were interpreting the neg/pos in the correct direction. 

We used a volcano plot to see that there is a lot more significant up regulation than down regulation. There is a point of data for each transcript. Our significant thresholds were set at p<0.05 and Log2FoldChange is >1. 

In the heat map. we can see how how large the differences of gene expression are for the top 20 most significantly differential expressed genes. We talked about how they were clustered by DevTemp but not clustered by FinalTemp. 

### 10/17/2024 Transcriptomics day 4 (contrasts)

Today we continuted working with the copopod dataset. We had to run the first two days of code to produce the objects we were working with today. We loaded the new package eulerr and most of today's work was in preparation to create a Eulerr plot. This is a plot with overlapping circles that looks like a venn diagram. It is used to show contrasts from an RNAseq experiment where each circle's radius increases as its number of differentially expressed genes (DEGs) increases. They overlap proportionally to the amount of DEGs shared.
We find the DEGs by reducing our contrast data frames to only include the rows that have a p value < 0.05.

### 10/22/2024 Transcriptomics day 5 
We created scatter plots of copepod contrasts allowing us to see how gene expression at a specific FinalTemp changed depending on the DevTemp (BASE vs 18). 
### 10/24/2024 Transcriptomics day 6
We began by finishing the scatter plots we started last class by adding in ablines. We then used grid.arrange function in the gridExtra package to combine the two plots into one figure. We used ggsave to export a png of the combined figure. The first step for this analysis is to import our data (counts matrix, conds data, and a new file which is our trait observation data)

We created a new script titled WGCNA. This script has the purpose of analyzing and visualizing gene correlation networks.

We ran a network topology analysis function then created plots of the power vs scale free topology and mean connectivity to pick a power that will create biologically relevant topology model. You want to pick an r^2 value > 0.8 without getting too close to zero on the power vs mean connectivity graph. 

The blockwise module commmand did not finish running and I should try again by increasing the amount of memory requested from the VACC to see if it runs. 
