# R script to use on hoffman
# module load R
# R 
# install R packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("gdsfmt")
#biocLite("SNPRelate")
# usage:
# module load R
# Rscript filtering_Step_4_ConvertToGDS.Hoffman.R [set the calldate and indir and infile manually, or can adjust script to make it a command line option.]
#load R packages
library(gdsfmt)
library(SNPRelate)

calldate=20180724 # date gt's were called in format YYYYMMDD (set this manually)
todaysdate=format(Sys.Date(),format="%Y%m%d")
SCRATCH="/u/flashscratch/a/ab08028"
indir=paste(SCRATCH,"/captures/vcf_filtered/",calldate,"/",sep="") # this is where your snp vcf file is and where you will save your gds file
infilePREFIX="snp_5_passingAllFilters_postMerge_raw_variants" # exclude the .vcf.gz suffix


#read vcf, and reformat to gds (this works with gzipped vcf file)

vcf.fn = paste(indir,infilePREFIX,".vcf.gz",sep="")

# took ~ 5 mins (only need to do once -- in future can just open the gds file)
snpgdsVCF2GDS(vcf.fn, paste(indir,infilePREFIX,".gds",sep=""), method="biallelic.only")

#summary -- write it out
# Open an output file
sink(file=paste(indir,infilePREFIX,".gds.summary.txt",sep=""))
snpgdsSummary(paste(indir,infilePREFIX,".gds",sep=""))
# close summary file
sink()

# you can then download the gds file and work with it in R on your home computer: for me download to: /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/snps_gds
