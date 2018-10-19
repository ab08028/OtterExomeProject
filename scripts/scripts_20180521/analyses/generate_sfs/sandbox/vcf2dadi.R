install.packages("devtools")
library(devtools)
install_github("thierrygosselin/radiator")
library("optparse")
library("radiator")
option_list = list(
  make_option(c("-vcf", "--vcf"), type="character", default=NULL, 
              help="path to vcf file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-pop", "--pop"), type="character", default=NULL, 
              help="population identifier", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

data=opt$vcf
pop=opt$vcf

vcf2dadi(data)