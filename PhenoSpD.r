install.packages("installr")
library(installr)
updateR()
source("https://bioconductor.org/biocLite.R")
biocLite("metaCCA")

#documentation
browseVignettes("metaCCA")
rm(list=ls(all=TRUE))

library(metaCCA)
data( package = 'metaCCA' )
print (N1)
View(S_XY_full_study1)
View(S_XX_study1)

S_YY_study1 = estimateSyy( S_XY = S_XY_full_study1 )

S_YY_study1 = estimateSyy( S_XY = S_XY_study1 )
S_YY_study1



data <- read.table("Google Drive/working space/Postdoc_year1_year2/projects/metaCCA/S_XY_full_shin_metabolites.txt",header=T,row.names=1)
#save(data, file = "Google Drive/working space/Postdoc_year1_year2/projects/metaCCA/S_XY_full_metabolites.RData")
#load()
#title <- read.table("Google Drive/working space/Postdoc_year1_year2/projects/metaCCA/metabolites.title.txt",header=F)
#data2 <- data[,-1]
#rownames(data2) <- data[,1]

#data frame to character
#title2[] <- lapply(title, as.character)
#names(data2) <- title2

#lower case to upper case
#install.packages("dplyr")
#library(dplyr)
#data2 <- mutate_each(data2, funs(toupper))

S_YY_shin_metabolites = estimateSyy( S_XY = data2 )

save(data2,file="Shin_metabolites_GWAS.Rda")
load("Shin_metabolites_GWAS.Rda")

install.packages("readxl")
library(readxl)
obs <- read_excel("~/Google Drive/working space/Postdoc_year1_year2/projects//metaCCA//phenotypic_correlation_comparison_obs_metaCCA.xlsx",1)
metaCCA <- read_excel("~/Google Drive/working space/Postdoc_year1_year2/projects//metaCCA//phenotypic_correlation_comparison_obs_metaCCA.xlsx",2)
obs2 <- obs[,-1]
rownames(obs2) <- obs[,1]

metaCCA2 <- metaCCA[,-1]
rownames(metaCCA2) <- metaCCA[,1]

#obs.list <- split(obs2, seq(nrow(obs2)))
#metaCCA.list <- split(metaCCA2, seq(nrow(metaCCA2)))

obs.list <- unlist(as.list(as.data.frame(t(obs2))))
metaCCA.list <- unlist(as.list(as.data.frame(t(metaCCA2))))

plot(obs.list,metaCCA.list,cex = 0.5)
abline(a=0, b=1,col="red")
abline(lm(metaCCA.list ~ obs.list),col="blue")

##############################SpD####################################
## To analyse your correlation matrix, execute this R script in the sample directory as your "correlation.matrix" file, using:
##within R:
## source("SpD.R")
## Results are written to 
## "matSpD.out"

source("SpD.R")


