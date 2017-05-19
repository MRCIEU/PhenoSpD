##THE SPD SECTION OF THE SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
##If you are estimating the number of independent traits, please cite: Nyholt DR (2004) A simple correction for multiple testing for SNPs in linkage disequilibrium with each other. Am J Hum Genet 74(4):765-769.

usePackage <- function(p) {
    if (!is.element(p, installed.packages()[,1]))
        install.packages(p, dep = TRUE)
    require(p, character.only = TRUE)
}

usePackage(optparse)

option_list = list(
  make_option(c("-s", "--phenocorr"), type="character", default=NULL, 
              help="file with the phenotypic correlation matrix", metavar="character"),

  make_option(c("-o", "--out"), type="character", default="pheno.corr.txt", 
              help="output of the multiple testing correction", metavar="character"),
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$phenocorr)){
  #print_help(opt_parser)
  stop("Please provide the file name of the phenotypic correlaiton matrix. \n", call.=FALSE)
}

##example commond line
#Rscript ./script/SpD.r --phenocorr example.pheno.corr.txt --out example.phenospd.txt


## Read in correlation matrix:    
# For multiple test correction the sign of the correlation is irrelevant (i.e., so we're best to input absolute values)
corr.matrix<-abs(read.table(opt$phenocorr))  

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

sink("opt$out")
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
  sink("opt$out")
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
