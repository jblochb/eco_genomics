getwd
setwd("/gpfs1/cl/pbio3990/Intro_to_R")
install.packages("tidyverse")
library(tidyverse)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiome")
Met<- read_csv("ThermalStressMetabolic.csv")
Met


