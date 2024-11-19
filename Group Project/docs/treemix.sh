#!/bin/bash

# cd command is used to call a directory
cd projects/eco_genomics/

# ll is used to view the files in a directory
ll

# -rw-r--r-- This is in front of each file and tells you the permissions
# the first character tells you if it is a directory (d) or file (-)
# Then the three character groups tell you the permission for owner, group, then everyone
# chmod o+x treemix_infile.txt this chmod command changes the permissions
# In this case, we are adding transformation access (x) to the owner (o)
# Another example:
# chmod o-x treemix_infile.txt 
# [jblochbe@vacc-user2 population_genomics]$ chmod u+x treemix_infile.txt
# Use the below command to zip a file 
gzip treemix_infile.txt
#The below command calls treemix from its directory then works on the input file and creates an output file
/gpfs1/cl/pbiosw/treemix-1.13/src/treemix -i treemix_infile.txt.gz -o treemix_output
# The output is a tree in Newick format