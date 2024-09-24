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

### 09/19/2024 - Diversity Differentiation

#Today we created the 02_Diversity_Differentiation_r script. We created a Manhattan plot showing the Fst value for each SNP. We were able to separate the chromosomes and see approximately where the centromere is (there are less SNPs there).

#This required reformatting our vcf file to the proper format and column titles that the qqman function requires.

#Some minor helpful things to remember: 
#If you want to change the type of a value in a dataframe (such as from character to numberic), you can set the column name equal to as.numberic(colname)

### 09/24/2024

#We made a histogram of Hs using ggplot and used the ggsave command to save it as a PDF
#We made a new R script for making a PCA. We also took notes on PCA. We had to filter the data to make sure that the SNPs were all separated by at least 500bp
