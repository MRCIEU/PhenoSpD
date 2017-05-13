###install.packages
#install.packages("installr")
#library(installr)
#updateR()


###estimate phenotypic correlation matrix using metaCCA###
#rm(list=ls(all=TRUE))
source("https://bioconductor.org/biocLite.R")
biocLite("metaCCA")
library(metaCCA)
data <- read.table("S_XY_full_shin_metabolites.txt",header=T,row.names=1)
S_YY_shin_metabolites = estimateSyy( S_XY = data)



##############################SpD####################################
## To analyse your correlation matrix, execute this R script in the sample directory as your "correlation.matrix" file, using:
##within R:
## source("SpD.R")
## Results are written to 
## "matSpD.out"

source("SpD.R")


