## Documenting Our Group Project Progress

### Thursday November 7th

We made out practice pitch to the class. We are going to work on piecing together the demographic history of the invasive Centaurea x monktonii complex in the United States. Some analyses we plan to use are treemix, PCA, approximate Bayesian computation, and STRUCTURE.

### Tuesday November 12th

We worked on organizing our ideas and figuring out what we want to work on. I chose to start with making a PCA of the capitula data, so we will be able to graph the PCs of this against the genomic data generated PCs. I was able to graph the capitula data using Steve's code. I will have to choose whether I want to use the 3 trait generated PCA for the future of the 9 trait generated one. \### Thursday November 14th Today I got the capitula PCs vs genomic PCs visualized by extracting the projections values and merging them into one data frame.

### Tuesday November 19th

I decided to create a data frame that has the data from PCs 1-3 for both the capitula generated data and the PC generated data in order to easily and quickly regenerate PCs based on the combination of any of these 6 data sets. I wrote this as a csv titled DF_for_PCA in my population_genomics folder. Treemix is my next task to work on, and it has now been downloaded onto the command line. Steve walked me through some basic bash coding so I was able to run the command and visualize the resulting tree. We used the vcf2others package to convert our filtered vcf file into a txt file suitable for treemix.

### Thursday November 21st

I kept working on treemix and was able to generate trees with varying amounts of migration events. I found it is easier to use the plot_tree() command within the treemix source code than to use a different R package to plot the Newick tree that treemix produces. We are thinking about the possibility of using an out group to establish a root of the tree so it can be interpreted photogenically.

### Tuesday December 3rd - December 5th 

Steve worked with me to merge our original vcf with the new six sample Centaurea solstitialis vcf that will serve as our out group. I then manually edited the metadata file I wrote as a csv using github to add the individual ids from the six new samples. When I pulled this from github back into RStudio, I was able to use it to merge the vcfs. I generated the out group trees and was surprised that there were many predicted gene flow events with the out group. For this reason, we decided to not use these trees for our conclusions. I then filtered much of the missingness out of the new C. solstitialis reference file, but the same error was occuring. 

### Finals week (December 9th - 13th)

We are writing our final paper as a group and making connections between our analyses. 