# Coding and data notes for the population genomics module

## Author: John Blochberger

### 09/10/2024 - Intro to Centaurea GBS data working with VCF files

Will be analyzing the GBS data from three regions (EU, NE, PNW) starting today with variant call format files (VCF's)

Accessed a fastq file on the command line by first navigating to the directory cd /gpfs1/cd/pbio3990 cd example_data/

to list the contents of the directory we use ls -l (can be abbreviated as ll) zcat file.gz \| head This is used to see the first few lines of data in a fastq file. To specify how many lines add (-n \#) after head

The \| symbol is the pipeline in bash coding language and it is used to send output of one function to another.

### 09/12/2024 - Viewing VCF files and talking about filtering 

spack samtools 
samtools tview (file)
space bar shifts right and delete key shifts left
white is a q mapping score of 30 and blue text has a q score of 20
  Use the question mark button to pop up a legend with what everything means

In bash coding, use cd .. to go back to previous directory