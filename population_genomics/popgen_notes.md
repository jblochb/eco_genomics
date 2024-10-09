# Coding and data notes for the population genomics module

## Author: John Blochberger

### 09/10/2024 - Intro to Centaurea GBS data working with VCF files

Will be analyzing the GBS data from three regions (EU, NE, PNW) starting today with variant call format files (VCF's)

Accessed a fastq file on the command line by first navigating to the directory cd /gpfs1/cd/pbio3990 cd example_data/

to list the contents of the directory we use ls -l (can be abbreviated as ll) zcat file.gz \| head This is used to see the first few lines of data in a fastq file. To specify how many lines add (-n \#) after head

The \| symbol is the pipeline in bash coding language and it is used to send output of one function to another.

### 09/12/2024 - Viewing VCF files and talking about filtering

spack load samtools samtools tview (file) space bar shifts right and delete key shifts left white is a q mapping score of 30 and blue text has a q score of 20 Use the question mark button to pop up a legend with what everything means

In bash coding, use cd .. to go back to previous directory

### 09/17/2024 - Filtering VCF files and performing diversity stats

There are two commands we can use to try to fix the error when graphing that sometimes occurs when coding on the vacc. They are options(bitmapType='cairo') and X11.options(type="cairo").

We learned about how genomic data is filtered. We filtered our vcf file by missingness, depth, and low frequency alleles.

-   It's important to have both a low and high depth filter. The low filter is used to remove SNPs that do not have enough coverage to produce an accurate base call and the high depth filter is used to remove SNPs that may have multiple DNA fragments from different genes mapping to them. This may be caused by gene paralogy.
-   The missingness filter is applied both at the SNP level and the individual level. This is used to remove SNPs that are not present in a significant amount of individuals and to remove individuals that do not have a significant amount of SNPs.
-   The low frequency allele filter is to remove alleles that are present in a very low percentage of the population such as 1%. These are removed because they have very information regarding phenotypes/genotypes at the population level.

### 09/19/2024 - Diversity Differentiation

Today we created the 02_Diversity_Differentiation_r script. We created a Manhattan plot showing the Fst value for each SNP. We were able to separate the chromosomes and see approximately where the centromere is (there are less SNPs there).

This required reformatting our vcf file to the proper format and column titles that the qqman function requires.

Some minor helpful things to remember: #If you want to change the type of a value in a dataframe (such as from character to numberic), you can set the column name equal to as.numberic(colname)

### 09/24/2024

We made a histogram of Hs using ggplot and used the ggsave command to save it as a PDF We made a new R script for making a PCA. We also took notes on PCA. We had to filter the data to make sure that the SNPs were all separated by at least 500bp to account for possible Linkage disequilibrium (LD).

We had to unzip the file outside of our working directory that we upload to github because it is so large. We then converted our vcf file to a geno file because this is the form the data should be in when running the pca command.

### 09/26/2024

We experimented with the settings for ggplotting a PC graph then saved it to our figures folder. I copied all the various edits we can make to the ggplot, so I can refer to this script in the future.

To find the amount of variance a PC axis accounts for, you can divide the eigen value of PC axis of interest by the sum of all the eigen values. The first axis will always account for the highest variance and 2nd the 2nd most variance, etc.

The plot that plots all the eigen values from greatest to least is known as a scree-plot.

Steve lectured on admixture mapping. He explained that the researcher chooses a K value which is the number of groups the program will separate the genomic data into. One individual can belong to multiple groups. This occurs when they have a history of introgression/admixture/hybridization. Calculate p, q, and 2pq. Iterate until you find the K value that produces the least error.

### 10/01/2024 - Admixture Analysis

Today we started with working on the 3 PCA Admixture R script from last class and created an admixture plot using the snmf function within the LEA package and imaging the analysis in a barplot. This function uses a .geno file. Within the snmf analysis, we used K = 1:10 to test a range of K values to see which might be the most accurate. We also set entropy = T to have the model use genomic data to test its strength. The K values we chose to investigate more were the ones in "the elbow of the graph." This means when the cross-entropy scores just begin to level off on the resulting plot. Interestingly, the plot of K values vs cross-entropy scores looks very similar when number of PCs is graphed against Eigen values. This is because we can think of our K value in terms of how many PCs we are using to categorize the data.

We also used the PCAdapt library for the first time. We filter out the minor frequency alleles with a frequency of less than or equal to 0.01. We also make sure to generate the PCA axes while accounting for LD. We also have to organize the SNPs in a data frame so they all have a column with an integer that designates the chromosome it is on. For the Manhattan command, all important data must be read as numeric and all NA values must be removed,

Other small notes

-   A Tibble is a tidyverse data frame.

-   You can label myK equal to any K value you are interested in. Using this throughout your code will allow you to investigate other K values faster.

-   dev.off() is required after exporting pdfs or finishing other commands such as par().


