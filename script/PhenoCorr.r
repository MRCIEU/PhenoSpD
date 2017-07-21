##############################################################################################################################
# PhenoSpD                                                                                                                   #
# Version 1.0.0                                                                                                              #
# (c) 2017 Jie Zheng and Tom Gaunt                                                                                           #
# University of Bristol, School of Social and Community Medicine                                                             #
# GNU General Public License v3                                                                                              #
##############################################################################################################################

usePackage <- function(p) {
    if (!is.element(p, installed.packages()[,1]))
        install.packages(p, dep = TRUE)
    require(p, character.only = TRUE)
}

usePackage(optparse)
usePackage(metaCCA)

option_list = list(
  make_option(c("-s", "--sumstats"), type="character", default=NULL, 
              help="file with mulitple GWAS results", metavar="character"),

  make_option(c("-o", "--out"), type="character", default="pheno.corr.txt", 
              help="output of the phenotypic correlation estimation", metavar="character"),


); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$sumstats)){
  #print_help(opt_parser)
  stop("Please provide the file name of the GWAS results. \n", call.=FALSE)
}

##example commond line
#Rscript ./script/PhenoCorr.r --input ./data/PhenoSpD_input_example.txt --out example.pheno.corr.txt

data <- read.table(opt$sumstats,header=T,row.names=1)

##estimate phenotypic correlation matrix using metaCCA 
S_YY = estimateSyy( S_XY = data)

##write the phenotypic correlation matrix into the output file
write.table(S_YY, file=opt$out, quote=F, sep="\t", row.names = F)

