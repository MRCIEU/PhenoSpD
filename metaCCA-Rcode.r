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



##############################matSpD####################################
## To analyse your correlation matrix, execute this R script in the sample directory as your "correlation.matrix" file, using:
## R CMD BATCH matSpD.R
## or from within R:
## source("matSpD.R")
##
## Results are written to 
## "matSpD.out"
## "matSpDunrotated.out"
## "matSpDvarimax.out"
## "matSpDpromax.out"
##
## Please cite
## Nyholt DR (2004) A simple correction for multiple testing for SNPs in linkage disequilibrium with each other. Am J Hum Genet 74(4):765-769.
## and
## http://neurogenetics.qimrberghofer.edu.au/matSpD/

## Read in correlation matrix:
#corr.matrix<-read.table("correlation.matrix")

#obs-shin
corr.matrix<-abs(read.table("/Users//Evelyn//Google Drive//working space//Postdoc_year1_year2//projects/metaCCA/corr.matrix.obs"))      # For multiple test correction the sign of the correlation is irrelevant (i.e., so we're best to input absolute values)

#metaCCA-shin
corr.matrix<-abs(read.table("/Users//Evelyn//Google Drive//working space//Postdoc_year1_year2//projects/metaCCA/corr.matrix.metaCCA"))      # For multiple test correction the sign of the correlation is irrelevant (i.e., so we're best to input absolute values)

#metaCCA-anthropometric
#corr.matrix<-abs(read.table("/Users//Evelyn//Google Drive//working space//Postdoc_year1_year2//projects/metaCCA/corr.matrix.anthropometric"))      # For multiple test correction the sign of the correlation is irrelevant (i.e., so we're best to input absolute values)

## Remove Duplicate Columns:
corr.matrix.RemoveDupCol <- corr.matrix[!duplicated((corr.matrix))]

## Remove Duplicate Rows:
corr.matrix.RemoveDupRow <- unique((corr.matrix.RemoveDupCol))

## Remove Redundant VAR Names:
VARnames.NonRedundant<-as.matrix(dimnames(corr.matrix.RemoveDupCol)[[2]])
colnames(VARnames.NonRedundant)<-"VAR"

evals<-eigen(t(corr.matrix.RemoveDupRow),symmetric=T)$values

oldV<-var(evals)
M<-length(evals)
L<-(M-1)
Meffold<-M*(1-(L*oldV/M^2))

if (evals == 1) { 
  oldV <- 0 
  Meffold <- M
}

labelevals<-array(dim=M)
for(col in 1:M) { labelevals[col]<-c(col) }
levals<-cbind(labelevals, evals)

newevals<-evals
for(i in 1:length(newevals)) { 
  if(newevals[i] < 0) { 
    newevals[i] <- 0
  }
}

newlevals<-cbind(labelevals, newevals)

newV<-var(newevals)
Meffnew<-M*(1-(L*newV/M^2))

if (evals == 1) { 
  newV <- 0 
  Meffnew <- M
}


##############################################################################################################################################

## Implement improved approach of Li and Ji. Heredity 2005 95:221-227

IntLinewevals<-newevals

for(i in 1:length(IntLinewevals)) {
  if(IntLinewevals[i] >= 1 ) {
    IntLinewevals[i] <- 1
  }
  if(IntLinewevals[i] < 1 ) {
    IntLinewevals[i] <- 0
  }
}

NonIntLinewevals <- newevals-floor(newevals)

MeffLi <- sum(NonIntLinewevals+IntLinewevals)

NewResultLitemp1<-c('Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005):')
NewResultLitemp2<-round(MeffLi,dig=4)
NewResultLi1<-matrix(NewResultLitemp1)
NewResultLi2<-matrix(NewResultLitemp2)
NewBonferroniLitemp<-c('Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%:',
                       ' ', 1-(0.95^(1/MeffLi)))
NewBonferroniLi<-matrix(NewBonferroniLitemp)
Separatortemp<-c(' ',
                 '--------------------------------------------------------------------------------',
                 ' ',                       ' ',
                 'USING THE MORE ACCURATE ESTIMATE OF THE Veff [VeffLi] PROPOSED BY LI AND JI (2005):')
Separator<-matrix(Separatortemp)
Messagetemp<-c(' ',
               'NB: I recommend using the Li and Ji (2005) estimate unless Veff < VeffLi. ')
Message<-matrix(Messagetemp)

##############################################################################################################################################



NewResulttemp<-c('Effective Number of Independent Variables [Veff] (-ve values set to zero):', round(Meffnew,dig=4))
NewBonferronitemp<-c('Experiment-wide Significance Threshold Required to Keep Type I Error Rate at 5%:', 1-(0.95^(1/Meffnew)))
NewEigenvaluestemp<-c('New Eigenvalues Associated with the Correlation Matrix:', round(newevals,dig=4))
NewVariancetemp<-c('Variance of the Eigenvalues (with -ve values set to zero):', round(newV,dig=4))

NewResult<-matrix(NewResulttemp)
NewBonferroni<-matrix(NewBonferronitemp)
NewEigenvalues<-matrix(NewEigenvaluestemp)
NewVariance<-matrix(NewVariancetemp)

Originaltemp<-c('Original (total) number of variables (V) after removing redundant (collinear) variables:',
                ' ', M)

OldEigenvalues1temp<-c(' ',
                       'For factor 1 to V, original eigenvalues associated with the correlation matrix:')
OldEigenvalues2temp<-round(newlevals,dig=4)

OldVariancetemp<-c(' ',
                   'Variance of the observed eigenvalues:', 
                   ' ', round(newV,dig=4))

OldResulttemp<-c(' ',
                 'Effective number of independent variables [Veff]:', 
                 ' ', round(Meffnew,dig=4))

OldBonferronitemp<-c(' ',
                     'Significance threshold required to keep Type I error rate at 5% (0.05/Veff):', 
                     ' ', 0.05/Meffnew)

Original<-matrix(Originaltemp)

OldResult<-matrix(OldResulttemp)
OldBonferroni<-matrix(OldBonferronitemp)
OldEigenvalues1<-matrix(OldEigenvalues1temp)
OldEigenvalues2<-OldEigenvalues2temp
OldVariance<-matrix(OldVariancetemp)

no.dimnames <- function(a) {
  ## Remove all dimension names from an array for compact printing.
  d <- list()
  l <- 0
  for(i in dim(a)) {
    d[[l <- l + 1]] <- rep("", i)
  }
  dimnames(a) <- d
  a
}

sink("/Users/Evelyn/Google Drive/working space/Postdoc_year1_year2//projects//metaCCA/matSpD.out")
print(no.dimnames(Original), quote=F)
print(no.dimnames(OldEigenvalues1), quote=F)
print(no.dimnames(OldEigenvalues2), quote=F)
print(no.dimnames(OldVariance), quote=F)
print(no.dimnames(OldResult), quote=F)
print(no.dimnames(OldBonferroni), quote=F)
print(no.dimnames(Separator), quote=F)
print(no.dimnames(NewResultLi1), quote=F)
print(no.dimnames(NewResultLi2), quote=F)
print(no.dimnames(NewBonferroniLi), quote=F)
print(no.dimnames(Message), quote=F)
sink()

Warningtemp<-c(' ',
               '### Warning ###: there were some negative eigenvalues!',
               'If the above results using negative eigenvalues are equivalent',
               'to the following - obtained by setting negative eigenvalues to zero -',
               'then the results should be fine.',
               ' ')
Warning<-matrix(Warningtemp)

if(any(evals < 0)) { 
  sink("/Users/Evelyn/Google Drive/working space/Postdoc_year1_year2//projects//metaCCA/matSpD.out")
  print(no.dimnames(Original), quote=F)
  print(no.dimnames(OldEigenvalues1), quote=F)
  print(no.dimnames(OldEigenvalues2), quote=F)
  print(no.dimnames(OldVariance), quote=F)
  print(no.dimnames(OldResult), quote=F)
  print(no.dimnames(OldBonferroni), quote=F)
  print(no.dimnames(Separator), quote=F)
  print(no.dimnames(NewResultLi1), quote=F)
  print(no.dimnames(NewResultLi2), quote=F)
  print(no.dimnames(NewBonferroniLi), quote=F)
  print(no.dimnames(Message), quote=F)
  sink()
}


#=============================================================================
# Perform Factor (Principal Components) Analysis on unrotated matrix

# Compute the singular-value decomposition
svdout <- La.svd(t(corr.matrix.RemoveDupRow))

# Get the eigenvalues from the singular-value decomposition
svdevals <- svdout$d

# Get the standard deviations of the principal components 
# (i.e., sqrt of the eigenvalues) from the singular-value decomposition
svdsdev <- sqrt(svdevals)

# Specify the Number of Factors to Extract
factors<-M

# Matrix of Variable Loadings (i.e., a matrix whose columns contain the eigenvectors)

Laloadings <- svdout$u

# Principal Component Coefficients for Unrotated Matrix

coefficients <- Laloadings [ ,1:factors ] %*% diag ( svdsdev[1:factors] )
urcoefficients <- coefficients

# Squared Principal Component Coefficients [eigenvectors] for Unrotated Matrix
urcoefficients2 <- urcoefficients^2

# Sum of Squared Principal Component Coefficients for Unrotated Matrix (i.e., eigenvalues)
urevals <- apply(urcoefficients, 2, function(x) sum(x^2) )

lurevals<-cbind(labelevals, urevals, urevals/cumsum(urevals)[M])

dimnames(lurevals)[[2]][1]<-"Eigenvalue"
dimnames(lurevals)[[2]][2]<-"Cumulative proportion of variance"

# Find the Maximum Value making up each Factor's Unrotated Eigenvalue
urmaxvalues<-array(dim=M)
for(col in 1:M) { urmaxvalues[col]<-c(round(max(urcoefficients2[,col]),dig=8))
}

# Designate which VAR Contributes the Most to each Unrotated Eigenvalue
urindependent<-matrix(data=0,nrow=M,ncol=M)
for(col in 1:M) {
  for(row in 1:M) {
    if(round(urcoefficients2[row,col],dig=8) == urmaxvalues[col])
      urindependent[row,col]<-1
    if(round(urcoefficients2[row,col],dig=8) == 0) urindependent[row,col]<-0
  }
}
#=============================================================================


OldurEigenvalues1temp<-c(' ',
                         '--------------------------------------------------------------------------------',
                         ' ',                       ' ',
                         'SELECT A SUBSET OF VARs WHILE OPTIMISING INFORMATION:',
                         ' ',
                         'For factor 1 to V, Eigenvalues and Proportion of Variance, with no Rotation:')
OldurEigenvalues2temp<-round(lurevals,dig=4)

Oldurcoefficients1temp<-c(' ',
                          'Principal component coefficients for unrotated matrix:',
                          ' - Columns represent factors (principal components) 1 to V',
                          ' - Rows represent Variable 1 to V',
                          ' ')
Oldurcoefficients2temp<-cbind(VARnames.NonRedundant,round(urcoefficients, dig=4))

dimnames(Oldurcoefficients2temp)[[2]][2:(M+1)]<-labelevals
dimnames(Oldurcoefficients2temp)[[1]]<-labelevals

OldurEigenLoadings1temp<-c(' ',
                           'Factor "loadings" with no rotation:',
                           ' - Columns represent factors 1 to V',
                           ' - Rows represent Variable 1 to V',
                           ' - Variables contributing the MOST to each unrotated factor are designated by a "1"',
                           ' ')


OldurEigenLoadings2temp<-cbind(VARnames.NonRedundant,urindependent)


dimnames(OldurEigenLoadings2temp)[[2]][2:(M+1)]<-labelevals
dimnames(OldurEigenLoadings2temp)[[1]]<-labelevals


OldurEigenLoadings3temp<-c(' => Select one VAR to represent either:',
                           '    i.   each factor,',
                           '    ii.  the factors with the largest Veff eigenvalues, or',
                           '    iii. the factors explaining a selected proportion of variance.')

OldurEigenvalues1<-matrix(OldurEigenvalues1temp)
OldurEigenvalues2<-OldurEigenvalues2temp
Oldurcoefficients1<-matrix(Oldurcoefficients1temp)
Oldurcoefficients2<-Oldurcoefficients2temp
OldurEigenLoadings1<-matrix(OldurEigenLoadings1temp)
OldurEigenLoadings2<-OldurEigenLoadings2temp
OldurEigenLoadings3<-matrix(OldurEigenLoadings3temp)

sink("/Users/Evelyn/Google Drive/working space/Postdoc_year1_year2//projects//metaCCA/matSpDunrotated.out")
print(no.dimnames(OldurEigenvalues1), quote=F)
print(no.dimnames(OldurEigenvalues2), quote=F)
print(no.dimnames(Oldurcoefficients1), quote=F)
print(Oldurcoefficients2, quote=F)
print(no.dimnames(OldurEigenLoadings1), quote=F)
print(OldurEigenLoadings2, quote=F)
print(no.dimnames(OldurEigenLoadings3), quote=F)
sink()


Warningtemp<-c(' ',
               '### Warning ###: there were some negative eigenvalues!',
               'If the above results using negative eigenvalues are equivalent',
               'to the following - obtained by setting negative eigenvalues to zero -',
               'then the results should be fine.',
               ' ')
Warning<-matrix(Warningtemp)

if(any(evals < 0)) { 
  sink("/Users/Evelyn/Google Drive/working space/Postdoc_year1_year2//projects//metaCCA/matSpDunrotated.out")
  print(no.dimnames(OldurEigenvalues1), quote=F)
  print(no.dimnames(OldurEigenvalues2), quote=F)
  print(no.dimnames(Oldurcoefficients1), quote=F)
  print(Oldurcoefficients2, quote=F)
  print(no.dimnames(OldurEigenLoadings1), quote=F)
  print(OldurEigenLoadings2, quote=F)
  print(no.dimnames(OldurEigenLoadings3), quote=F)
  
  sink()
}


#=============================================================================
# Perform VARIMAX Rotation to Sharpen VAR Loadings to Particular Eigenvectors

# Compute the singular-value decomposition
svdout <- La.svd(t(corr.matrix.RemoveDupRow))

# Get the eigenvalues from the singular-value decomposition
svdevals <- svdout$d

# Get the standard deviations of the principal components 
# (i.e., sqrt of the eigenvalues) from the singular-value decomposition
svdsdev <- sqrt(svdevals)

# Specify the Number of Factors to Extract
factors<-M

# Matrix of Variable Loadings (i.e., a matrix whose columns contain the eigenvectors)

Laloadings <- svdout$u

# Principal Component Coefficients for Unrotated Matrix

coefficients <- Laloadings [ ,1:factors ] %*% diag ( svdsdev[1:factors] )

# Principal Component Coefficients for Varimax Rotated Matrix
vrcoefficients <- varimax ( coefficients ) $loadings 

# Squared Principal Component Coefficients [eigenvectors] for Varimax Rotated Matrix
vrcoefficients2 <- vrcoefficients^2

# Sum of Squared Principal Component Coefficients for Varimax Rotated Matrix (i.e., eigenvalues)
vrevals <- apply(vrcoefficients, 2, function(x) sum(x^2) )

lvrevals<-cbind(labelevals, vrevals, vrevals/cumsum(vrevals)[M])

dimnames(lvrevals)[[2]][1]<-"Eigenvalue"
dimnames(lvrevals)[[2]][2]<-"Cumulative proportion of variance"

# Find the Maximum Value making up each Factor's Rotated Eigenvalue
vrmaxvalues<-array(dim=M)
for(col in 1:M) { vrmaxvalues[col]<-c(round(max(vrcoefficients2[,col]),dig=8))
}

# Designate which VAR Contributes the Most to each Rotated Eigenvalue
vrindependent<-matrix(data=0,nrow=M,ncol=M)
for(col in 1:M) {
  for(row in 1:M) {
    if(round(vrcoefficients2[row,col],dig=8) == vrmaxvalues[col])
      vrindependent[row,col]<-1
    if(round(vrcoefficients2[row,col],dig=8) == 0) vrindependent[row,col]<-0
  }
}
#=============================================================================


OldvrEigenvalues1temp<-c(' ',
                         '--------------------------------------------------------------------------------',
                         ' ',                       ' ',
                         'SELECT A SUBSET OF VARs WHILE OPTIMISING INFORMATION:',
                         ' ',
                         'For factor 1 to V, Eigenvalues and Proportion of Variance, after Varimax Rotation:')
OldvrEigenvalues2temp<-round(lvrevals,dig=4)

Oldvrcoefficients1temp<-c(' ',
                          'Principal component coefficients for varimax-rotated matrix:',
                          ' - Columns represent factors (principal components) 1 to V',
                          ' - Rows represent Variable 1 to V',
                          ' ')
Oldvrcoefficients2temp<-cbind(VARnames.NonRedundant,round(vrcoefficients, dig=4))

dimnames(Oldvrcoefficients2temp)[[2]][2:(M+1)]<-labelevals
dimnames(Oldvrcoefficients2temp)[[1]]<-labelevals

OldvrEigenLoadings1temp<-c(' ',
                           'Factor "loadings" after varimax rotation:',
                           ' - Columns represent factors 1 to V',
                           ' - Rows represent Variable 1 to V',
                           ' - Variables contributing the MOST to each rotated factor are designated by a "1"',
                           ' ')


OldvrEigenLoadings2temp<-cbind(VARnames.NonRedundant,vrindependent)


dimnames(OldvrEigenLoadings2temp)[[2]][2:(M+1)]<-labelevals
dimnames(OldvrEigenLoadings2temp)[[1]]<-labelevals


OldvrEigenLoadings3temp<-c(' => Select one VAR to represent either:',
                           '    i.   each factor,',
                           '    ii.  the factors with the largest Veff eigenvalues, or',
                           '    iii. the factors explaining a selected proportion of variance.')

OldvrEigenvalues1<-matrix(OldvrEigenvalues1temp)
OldvrEigenvalues2<-OldvrEigenvalues2temp
Oldvrcoefficients1<-matrix(Oldvrcoefficients1temp)
Oldvrcoefficients2<-Oldvrcoefficients2temp
OldvrEigenLoadings1<-matrix(OldvrEigenLoadings1temp)
OldvrEigenLoadings2<-OldvrEigenLoadings2temp
OldvrEigenLoadings3<-matrix(OldvrEigenLoadings3temp)

sink("/Users/Evelyn/Google Drive/working space/Postdoc_year1_year2//projects//metaCCA/matSpDvarimax.out")
print(no.dimnames(OldvrEigenvalues1), quote=F)
print(no.dimnames(OldvrEigenvalues2), quote=F)
print(no.dimnames(Oldvrcoefficients1), quote=F)
print(Oldvrcoefficients2, quote=F)
print(no.dimnames(OldvrEigenLoadings1), quote=F)
print(OldvrEigenLoadings2, quote=F)
print(no.dimnames(OldvrEigenLoadings3), quote=F)
sink()


Warningtemp<-c(' ',
               '### Warning ###: there were some negative eigenvalues!',
               'If the above results using negative eigenvalues are equivalent',
               'to the following - obtained by setting negative eigenvalues to zero -',
               'then the results should be fine.',
               ' ')
Warning<-matrix(Warningtemp)

if(any(evals < 0)) { 
  sink("/Users/Evelyn/Google Drive/working space/Postdoc_year1_year2//projects//metaCCA/matSpDvarimax.out")
  print(no.dimnames(OldvrEigenvalues1), quote=F)
  print(no.dimnames(OldvrEigenvalues2), quote=F)
  print(no.dimnames(Oldvrcoefficients1), quote=F)
  print(Oldvrcoefficients2, quote=F)
  print(no.dimnames(OldvrEigenLoadings1), quote=F)
  print(OldvrEigenLoadings2, quote=F)
  print(no.dimnames(OldvrEigenLoadings3), quote=F)
  
  sink()
}



#=============================================================================
# Perform PROMAX Rotation to Sharpen VAR Loadings to Particular Eigenvectors

# Compute the singular-value decomposition
svdout <- La.svd(t(corr.matrix.RemoveDupRow))

# Get the eigenvalues from the singular-value decomposition
svdevals <- svdout$d

# Get the standard deviations of the principal components 
# (i.e., sqrt of the eigenvalues) from the singular-value decomposition
svdsdev <- sqrt(svdevals)

# Specify the Number of Factors to Extract
factors<-M

# Matrix of Variable Loadings (i.e., a matrix whose columns contain the eigenvectors)

Laloadings <- svdout$u

# Principal Component Coefficients for Unrotated Matrix

coefficients <- Laloadings [ ,1:factors ] %*% diag ( svdsdev[1:factors] )

# Principal Component Coefficients for Promax Rotated Matrix
prcoefficients <- promax ( coefficients ) $loadings 

# Squared Principal Component Coefficients [eigenvectors] for Promax Rotated Matrix
prcoefficients2 <- prcoefficients^2

# Sum of Squared Principal Component Coefficients for Promax Rotated Matrix (i.e., eigenvalues)
prevals <- apply(prcoefficients, 2, function(x) sum(x^2) )

lprevals<-cbind(labelevals, prevals, prevals/cumsum(prevals)[M])

dimnames(lprevals)[[2]][1]<-"Eigenvalue"
dimnames(lprevals)[[2]][2]<-"Cumulative proportion of variance"

# Find the Maximum Value making up each Factor's Rotated Eigenvalue
prmaxvalues<-array(dim=M)
for(col in 1:M) { prmaxvalues[col]<-c(round(max(prcoefficients2[,col]),dig=8))
}

# Designate which VAR Contributes the Most to each Rotated Eigenvalue
prindependent<-matrix(data=0,nrow=M,ncol=M)
for(col in 1:M) {
  for(row in 1:M) {
    if(round(prcoefficients2[row,col],dig=8) == prmaxvalues[col])
      prindependent[row,col]<-1
    if(round(prcoefficients2[row,col],dig=8) == 0) prindependent[row,col]<-0
  }
}
#=============================================================================


OldprEigenvalues1temp<-c(' ',
                         '--------------------------------------------------------------------------------',
                         ' ',                       ' ',
                         'SELECT A SUBSET OF VARs WHILE OPTIMISING INFORMATION:',
                         ' ',
                         'For factor 1 to V, Eigenvalues and Proportion of Variance, after Promax Rotation:')
OldprEigenvalues2temp<-round(lprevals,dig=4)

Oldprcoefficients1temp<-c(' ',
                          'Principal component coefficients for promax-rotated matrix:',
                          ' - Columns represent factors (principal components) 1 to V',
                          ' - Rows represent Variable 1 to V',
                          ' ')
Oldprcoefficients2temp<-cbind(VARnames.NonRedundant,round(prcoefficients, dig=4))

dimnames(Oldprcoefficients2temp)[[2]][2:(M+1)]<-labelevals
dimnames(Oldprcoefficients2temp)[[1]]<-labelevals

OldprEigenLoadings1temp<-c(' ',
                           'Factor "loadings" after promax rotation:',
                           ' - Columns represent factors 1 to V',
                           ' - Rows represent Variable 1 to V',
                           ' - Variables contributing the MOST to each rotated factor are designated by a "1"',
                           ' ')


OldprEigenLoadings2temp<-cbind(VARnames.NonRedundant,prindependent)


dimnames(OldprEigenLoadings2temp)[[2]][2:(M+1)]<-labelevals
dimnames(OldprEigenLoadings2temp)[[1]]<-labelevals


OldprEigenLoadings3temp<-c(' => Select one VAR to represent either:',
                           '    i.   each factor,',
                           '    ii.  the factors with the largest Veff eigenvalues, or',
                           '    iii. the factors explaining a selected proportion of variance.')

OldprEigenvalues1<-matrix(OldprEigenvalues1temp)
OldprEigenvalues2<-OldprEigenvalues2temp
Oldprcoefficients1<-matrix(Oldprcoefficients1temp)
Oldprcoefficients2<-Oldprcoefficients2temp
OldprEigenLoadings1<-matrix(OldprEigenLoadings1temp)
OldprEigenLoadings2<-OldprEigenLoadings2temp
OldprEigenLoadings3<-matrix(OldprEigenLoadings3temp)

sink("/Users/Evelyn/Google Drive/working space/Postdoc_year1_year2//projects//metaCCA/matSpDpromax.out")
print(no.dimnames(OldprEigenvalues1), quote=F)
print(no.dimnames(OldprEigenvalues2), quote=F)
print(no.dimnames(Oldprcoefficients1), quote=F)
print(Oldprcoefficients2, quote=F)
print(no.dimnames(OldprEigenLoadings1), quote=F)
print(OldprEigenLoadings2, quote=F)
print(no.dimnames(OldprEigenLoadings3), quote=F)
sink()


Warningtemp<-c(' ',
               '### Warning ###: there were some negative eigenvalues!',
               'If the above results using negative eigenvalues are equivalent',
               'to the following - obtained by setting negative eigenvalues to zero -',
               'then the results should be fine.',
               ' ')
Warning<-matrix(Warningtemp)

if(any(evals < 0)) { 
  sink("/Users/Evelyn/Google Drive/working space/Postdoc_year1_year2//projects//metaCCA/matSpDpromax.out")
  print(no.dimnames(OldprEigenvalues1), quote=F)
  print(no.dimnames(OldprEigenvalues2), quote=F)
  print(no.dimnames(Oldprcoefficients1), quote=F)
  print(Oldprcoefficients2, quote=F)
  print(no.dimnames(OldprEigenLoadings1), quote=F)
  print(OldprEigenLoadings2, quote=F)
  print(no.dimnames(OldprEigenLoadings3), quote=F)
  
  sink()
}
