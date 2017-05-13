###estimate phenotypic correlation matrix using metaCCA###
#rm(list=ls(all=TRUE))
source("https://bioconductor.org/biocLite.R")
biocLite("metaCCA")
library(metaCCA)
data <- read.table("S_XY_full_shin_metabolites.txt",header=T,row.names=1)
S_YY_shin_metabolites = estimateSyy( S_XY = data)
