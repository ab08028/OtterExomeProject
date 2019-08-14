############# In the previous step (7b-i), the script pulled out the scaffolds of the ferret genome for which you have cds information in your input file. It then binned those scaffolds into bins of a specified bin size.  For the specified individual, it got totals of differnet categories of site per bin. Many of these bins will be empty because they don't contain cds sites, so they will be dropped. Each remaining bin will have the GPs from ANGSD summed up (for sites *passing depth and GP filters*) to get the sums of HomRef, homAlt and het sites. It did this for all sites and then for transversions only. It wrote out these bin totals which are the input for this step. In this step (7b-ii), it will resample (with replacement) that individual's bins until a threshold of sites is reached (which you specify, should be average callable cds sites across all your individuals).  It will do this a specified number of times (numBoot) to generate bootstrap distributions for each individual. 

# This is an improvement over my old way of doing it that just resampled sites. By resampling over bins of the genome we are getting at LD structure better. 

# Note this sums up combo annotations because it assumes --pick has been used! (so if it's missense,splice_something) it still gets to be missense; synonymous,splice is also syn.
################## left to do: somewhere need to check for uniqueness (maybe before the script) #### 
#transversions=c('A,C','C,A','A,T','T,A','C,G','G,C','G,T','T,G')
#BiocManager::install("plyranges")
require(GenomicRanges)
require(dplyr)
require(plyranges)
#require(bootstrap)
# install plyranges (dplyr to work with Granges)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("plyranges")
## okay this works pretty well; want to parallelize and try running on hoffman and be able to feed in files; and make sure point estimates are inside there somewhere.
#require(plyranges) # this is needed so you can use group_by with GRanges objects!!
#minDepth=1 # make this match whatever I used to get point estimate
#minGP=0.95 # make this match whatever I used to get point estimate
#binsize=100000 # start with 100kb
#numBoots=5 # do boots a number of times
#SitesToDraw=30000 # <-- make this the avg callable cds sites across all individuals
#out.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/bootStrapRegions/"
library("optparse")
option_list = list(
  make_option(c("--infile"), type="character", default=NULL, 
              help="path to your input file result of step 7b-i which generates per indivdiual bins across the genome with the GPs of each category (syn, mis, sg) summed up. This is the df indAllBins from previous step", metavar="file"),
  make_option(c("--outdir"), type="character", default=NULL, 
              help="path to outdir", metavar="dir"),
  make_option(c("--outPREFIX"), type="character", default=NULL, 
              help="outfilePrefix", metavar="prefix"),
  make_option(c("--numBoots"), type="numeric", default=NULL, 
              help="Number of bootstraps to perform per individual", metavar="numeric"),
  make_option(c("--avgSitesToDraw"), type="numeric", default=NULL, 
              help="Number of sites to draw per individual (will be approximately this many, not exactly) due to varying bin size. Based on average 'callable' cds sites calculated elsewhere.", metavar="numeric"),
  make_option(c("--indNum"), type="numeric", default=NULL, 
              help="Individual number assigned by ANGSD, starts at 0 (for my study it's 0-8))", metavar="numeric")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

infile=opt$infile # output by step a-i - paste(out.dir,"/",outPREFIX,"Ind.",ind,".sumsPerBin.txt",sep="")
out.dir=opt$outdir
outPREFIX=opt$outPREFIX
numBoots=as.numeric(opt$numBoots) # do boots a number of times
SitesToDraw=as.numeric(opt$avgSitesToDraw)
ind=as.numeric(opt$indNum)
# read in results of step 7b-i:
indAllBins = read.table(infile,header=T,sep="\t")

# test:
#outPREFIX="TEST"
#numBoots=5
#SitesToDraw=30000
#infile="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/bootStrapRegions/TESTInd.1.sumsPerBin.txt"
#out.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/bootStrapRegions/"
##### randomly sample bins from an individual until some number is reached
allBoots=data.frame()
print("starting resampling bootstraps")
for(i in seq(1,numBoots)){
  print(c("starting bootstrap: ",i))
  runningTotal=0
  numberOfDraws=0
  binNums=unique(indAllBins$binNum)
  ###### if this has been through "pick" then can use grepl because 
  
  # can add other categories if I want
  MISdf=data.frame(ind=ind,bootNum=i,homRef=0,het=0,homAlt=0,category="missense",siteType="Ti+Tv")
  SYNdf=data.frame(ind=ind,bootNum=i,homRef=0,het=0,homAlt=0,category="synonymous",siteType="Ti+Tv")
  SGdf=data.frame(ind=ind,bootNum=i,homRef=0,het=0,homAlt=0,category="stop_gained",siteType="Ti+Tv")
  # tv only:
  MISdfTV=data.frame(ind=ind,bootNum=i,homRef=0,het=0,homAlt=0,category="missense",siteType="TvOnly")
  SYNdfTV=data.frame(ind=ind,bootNum=i,homRef=0,het=0,homAlt=0,category="synonymous",siteType="TvOnly")
  SGdfTV=data.frame(ind=ind,bootNum=i,homRef=0,het=0,homAlt=0,category="stop_gained",siteType="TvOnly")
  while(runningTotal < SitesToDraw){
    randomBin=sample(binNums,1) # draw a bin number -- will have replacement because you're not removing the number
    bin=indAllBins[indAllBins$binNum==randomBin,]
    numberOfDraws=numberOfDraws+1
    runningTotal = runningTotal+unique(bin$totalCallableSitesPerBin)
    if((runningTotal+unique(bin$totalCallableSitesPerBin)) < SitesToDraw){
      MISdf$homRef= MISdf$homRef + sum(bin[grepl("missense_variant",bin$Consequence) & bin$sites=="Ti+Tv",]$sumHomRef)
      MISdf$het= MISdf$het +sum(bin[grepl("missense_variant",bin$Consequence) & bin$sites=="Ti+Tv",]$sumHet)
      MISdf$homAlt= MISdf$homAlt +sum(bin[grepl("missense_variant",bin$Consequence) & bin$sites=="Ti+Tv",]$sumHomAlt)
      # synonymous:
      SYNdf$homRef= SYNdf$homRef + sum(bin[grepl("synonymous_variant",bin$Consequence) & bin$sites=="Ti+Tv",]$sumHomRef)
      SYNdf$het= SYNdf$het +sum(bin[grepl("synonymous_variant",bin$Consequence) & bin$sites=="Ti+Tv",]$sumHet)
      SYNdf$homAlt= SYNdf$homAlt +sum(bin[grepl("synonymous_variant",bin$Consequence) & bin$sites=="Ti+Tv",]$sumHomAlt)
      # stop_gained:
      SGdfTV$homRef= SGdfTV$homRef + sum(bin[grepl("stop_gained",bin$Consequence) & bin$sites=="Ti+Tv",]$sumHomRef)
      SGdfTV$het= SGdfTV$het +sum(bin[grepl("stop_gained",bin$Consequence) & bin$sites=="Ti+Tv",]$sumHet)
      SGdfTV$homAlt= SGdfTV$homAlt +sum(bin[grepl("stop_gained",bin$Consequence) & bin$sites=="Ti+Tv",]$sumHomAlt)
      # transversions:
      MISdfTV$homRef= MISdfTV$homRef + sum(bin[grepl("missense_variant",bin$Consequence) & bin$sites=="TvOnly",]$sumHomRef)
      MISdfTV$het= MISdfTV$het +sum(bin[grepl("missense_variant",bin$Consequence) & bin$sites=="TvOnly",]$sumHet)
      MISdfTV$homAlt= MISdfTV$homAlt +sum(bin[grepl("missense_variant",bin$Consequence) & bin$sites=="TvOnly",]$sumHomAlt)
      # synonymous:
      SYNdfTV$homRef= SYNdfTV$homRef + sum(bin[grepl("synonymous_variant",bin$Consequence) & bin$sites=="TvOnly",]$sumHomRef)
      SYNdfTV$het= SYNdfTV$het +sum(bin[grepl("synonymous_variant",bin$Consequence) & bin$sites=="TvOnly",]$sumHet)
      SYNdfTV$homAlt= SYNdfTV$homAlt +sum(bin[grepl("synonymous_variant",bin$Consequence) & bin$sites=="TvOnly",]$sumHomAlt)
      # stop_gained:
      SGdfTV$homRef= SGdfTV$homRef + sum(bin[grepl("stop_gained",bin$Consequence) & bin$sites=="TvOnly",]$sumHomRef)
      SGdfTV$het= SGdfTV$het +sum(bin[grepl("stop_gained",bin$Consequence) & bin$sites=="TvOnly",]$sumHet)
      SGdfTV$homAlt= SGdfTV$homAlt +sum(bin[grepl("stop_gained",bin$Consequence) & bin$sites=="TvOnly",]$sumHomAlt)
      
    }}
  # update with final running totals
  MISdf$totalCDSSitesInBoot <- runningTotal
  SYNdf$totalCDSSitesInBoot <- runningTotal
  SGdf$totalCDSSitesInBoot <- runningTotal
  MISdfTV$totalCDSSitesInBoot <- runningTotal
  SYNdfTV$totalCDSSitesInBoot <- runningTotal
  SGdfTV$totalCDSSitesInBoot <- runningTotal
  
  allBoots=rbind(allBoots,MISdf,MISdfTV,SYNdf,SYNdfTV,SGdf,SGdfTV) # update allBoots
}
write.table(allBoots,paste(out.dir,"/",outPREFIX,".Ind.",ind,".allBoots.txt",sep=""),col.names = T,row.names = F,quote=F)

# not writing out all inds all boots -- keeping separate per ind to parallelize