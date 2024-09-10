# Coding and data notes for the population genomics module

## Author: John Blochberger

### 09/10/2024 - Intro to Centaurea GBS data working with VCF files

Will be analyzing the GBS data from three regions (EU, NE, PNW) starting today with variant call format files (VCF's)

Accessed a fastq file on the command line by first navigating to the directory cd /gpfs1/cd/pbio3990 cd example_data/

to list the contents of the directory we use ls -l (can be abbreviated as ll)
zcat file.gz | head 
    This is used to see the first few lines of data in a fastq file. To specify     how many lines add (-n \#) after head

The | symbol is the pipeline in bash coding language and it is used to send output of one function to another.
