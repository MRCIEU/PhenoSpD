##THE SPD SECTION OF THE SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS OF SNPSpD AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
##If you are estimating the number of independent traits, please cite: Nyholt DR (2004) A simple correction for multiple testing for SNPs in linkage disequilibrium with each other. Am J Hum Genet 74(4):765-769.

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
    
  make_option(c("-s", "--phenocorr"), type="character", default=NULL, 
              help="phenotypic correlation matrix of multiple traits", metavar="character"),
    
  make_option(c("-o", "--out"), type="character", default="phenospd", 
              help="output of the phenotypic correlation estimation", metavar="character"),


); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$sumstats)==FALSE & is.null(opt$sumstats)==FALSE){
  print_help(opt_parser)
  stop("Please provide either the GWAS results file or the phenotypic correlaiton matrix file (you are providing both now). \n", call.=FALSE)
}

if (is.null(opt$sumstats)==TRUE & is.null(opt$sumstats)==TRUE){
  print_help(opt_parser)
  stop("No input file provided. \n", call.=FALSE)
}


##example commond line
#Rscript ./script/phenospd.r --sumstats ./data/PhenoSpD_input_example.txt --out example
#Rscript ./script/phenospd.r --phenocorr ./data/LD-Hub_phenotypic_correlation_221x221.txt --out example

if (is.null(opt$sumstats)==FALSE & is.null(opt$sumstats)==TRUE){
  data <- read.table(opt$sumstats,header=T,row.names=1)
  ##estimate phenotypic correlation matrix using metaCCA 
  S_YY = estimateSyy( S_XY = data)
  ##write the phenotypic correlation matrix into the output file
  out_phenocorr <-paste0(opt$out,"_corr.txt")
  write.table(S_YY, file=out_phenocorr, quote=F, sep="\t", row.names = F)
  ##estimate number of independent traits and write the results into the output file
  source(./script/SpD.r)
  out_ntests <-paste0(opt$out,"_ntests.txt")
  SpD(phenocorr=S_YY, out=out_ntests) 
}

if (is.null(opt$sumstats)==TRUE & is.null(opt$sumstats)==FALSE){
  ##estimate number of independent traits and write the results into the output file
  source(./script/SpD.r) 
  out_ntests <-paste0(opt$out,"_ntests.txt")
  SpD(phenocorr=S_YY, out=out_ntests)      
  #system(paste0("Rscript ./script/phenospd.r --phenocorr ",out_phenocorr," --out ",opt$out,"_ntests.txt")) 
}



