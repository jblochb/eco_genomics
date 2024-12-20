library(ggfortify)
library(tidyverse)
options(bitmapType = "cairo")
X11.options(type="cairo")
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/traits/capitula/")
setwd("~/")
dat <- read.csv("PNW_EU_NE_capitulummeasurements.csv", header=T)

str(dat)
head(dat)
dat = dat %>%
  mutate(id = paste(Pop, IndID, sep="_"), .before=5) %>%
  group_by(id) %>%
  slice_head(n=1)

head(dat)

summary(dat)
?prcomp
?stats

pca1 <- prcomp(dat[,c(7:20)], center=T, scale=T)
View(pca1)
pca1.df <- as.data.frame(pca1$x)
pca2 <- prcomp(dat[,c(16,18,19)], center=T, scale=T)
?prcomp
?autoplot
autoplot(pca1, data = dat, colour = 'Region',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

autoplot(pca2, data = dat, colour = 'Region',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)


citation()
