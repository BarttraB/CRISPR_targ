
###install bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite()
require("Biostrings")

library("Biostrings")

setwd("/Volumes/Storage/R/CRISPR_TargFind")

s <- readDNAStringSet("E.coli-MG1655-genome.fasta")

s1 <- DNAStringSet(s, start=2, end=20)
### Hamming distance
h=sum(x1 != x2);

w=t(a)%*%a
