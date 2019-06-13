# R script to use on hoffman
# module load R
# R 
# install R packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")
# usage:
# module load R
# Rscript filtering_Step_4_ConvertToGDS.Hoffman.R
# input: plink bed (binary ped) file; NOTE this is not a UCSC coordidate bed file, no relation 
# can be run in the shell (will take ~ 10 min)
#load R packages
require(gdsfmt)
require(SNPRelate)
require("optparse")
option_list = list(
  make_option(c("-PlinkBed", "--PlinkBed"), type="character", default=NULL, 
              help="path to bed file from plink", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="path to outputdir [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

outdir=opt$outdir
bed.fn=opt$PlinkBed
# note that basename will get you just the file name from a path ; dirname will get you just the dir
bed.fn.prefix=strsplit(basename(bed.fn),".bed") # strip off end 
outfile=paste(outdir,"/",bed.fn.prefix,".gds",sep="")

# took ~ 5 mins (only need to do once -- in future can just open the gds file)
snpgdsVCF2GDS(bed.fn, outfile, method="biallelic.only")

#summary -- write it out
# Open an output file
sink(file=paste(outfile,".summary.txt",sep=""))
snpgdsSummary(outfile)
# close summary file
sink()

# you can then download the gds file and work with it in R on your home computer: for me download to: /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/snps_gds
