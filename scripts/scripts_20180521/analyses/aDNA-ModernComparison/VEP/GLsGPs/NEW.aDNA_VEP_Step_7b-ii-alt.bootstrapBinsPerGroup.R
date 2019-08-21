############# 
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
              help="Number of sites to draw per individual (will be approximately this many, not exactly) due to varying bin size. Based on average 'callable' cds sites calculated elsewhere.", metavar="numeric")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

infile=opt$infile # output by step a-i - paste(out.dir,"/",outPREFIX,"Ind.",ind,".sumsPerBin.txt",sep="")
out.dir=opt$outdir
outPREFIX=opt$outPREFIX
numBoots=as.numeric(opt$numBoots) # do boots a number of times
SitesToDraw=as.numeric(opt$avgSitesToDraw)
# read in results of step 7b-i:
allBins = read.table(infile,header=T,sep="\t")
######### for testing #####
#allBins <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/bootStrapRegions/sandbox/angsdOut.mappedTomfur.Bins.GPs.ProbCutoff.0.95.DepthCutoff.2.minInd.1.20190701-highcov-AFprior-MajorMinor4.Modern.Ancient.AvgsPerGroup.PerBin.txt",header=T)
#head(allBins)
#numBoots=5
SitesToDraw=30000

##### randomly sample bins from an individual until some number is reached
allBoots=data.frame()
print("starting resampling bootstraps")
# do ancient and modern distributions separately
for(group in c("Ancient","Modern")){
  # select out ancient or modern rows: 
  GroupBoots <- data.frame()
  GroupBins <- allBins[allBins$group==group,]
  for(i in seq(1,numBoots)){
    print(c("starting bootstrap: ",i))
    runningTotal=0
    numberOfDraws=0
    binNums=unique(GroupBins$binNum)
    ###### if this has been through "pick" then can use grepl because 
    
    # can add other categories if I want
    MISdf=data.frame(group=group,bootNum=i,homRef=as.numeric(0),het=as.numeric(0),homAlt=as.numeric(0),category="missense",siteType="Ti+Tv",stringsAsFactors = F)
    SYNdf=data.frame(group=group,bootNum=i,homRef=as.numeric(0),het=as.numeric(0),homAlt=as.numeric(0),category="synonymous",siteType="Ti+Tv")
    SGdf=data.frame(group=group,bootNum=i,homRef=as.numeric(0),het=as.numeric(0),homAlt=as.numeric(0),category="stop_gained",siteType="Ti+Tv")
    # tv only:
    MISdfTV=data.frame(group=group,bootNum=i,homRef=as.numeric(0),het=as.numeric(0),homAlt=as.numeric(0),category="missense",siteType="TvOnly")
    SYNdfTV=data.frame(group=group,bootNum=i,homRef=as.numeric(0),het=as.numeric(0),homAlt=as.numeric(0),category="synonymous",siteType="TvOnly")
    SGdfTV=data.frame(group=group,bootNum=i,homRef=as.numeric(0),het=as.numeric(0),homAlt=as.numeric(0),category="stop_gained",siteType="TvOnly")
    while(runningTotal < SitesToDraw){
      randomBin=sample(binNums,1) # draw a bin number -- will have replacement because you're not removing the number
      bin=GroupBins[GroupBins$binNum==randomBin,]
      numberOfDraws=numberOfDraws+1
      runningTotal = runningTotal+unique(bin$averageCalledSitesPerBin)
      if((runningTotal+unique(bin$averageCalledSitesPerBin)) < SitesToDraw){
        MISdf$homRef= MISdf$homRef + sum(bin[grepl("missense_variant",bin$Consequence) & bin$sites=="Ti+Tv",]$avgRescaledHomRef)
        MISdf$het= MISdf$het +sum(bin[grepl("missense_variant",bin$Consequence) & bin$sites=="Ti+Tv",]$avgRescaledHet)
        MISdf$homAlt= MISdf$homAlt +sum(bin[grepl("missense_variant",bin$Consequence) & bin$sites=="Ti+Tv",]$avgRescaledHomAlt)
        # synonymous: ** note syn treated diff than mis or sg here to match with filter vep!** to match with point estimates from filter_vep want to make sure it's a) canonical and b) only equals synonymous_variant not splice_region_variant,synonymous_variant (this is a choice -- doesn't really matter, just want to be consistent.)
        SYNdf$homRef= SYNdf$homRef + sum(bin[bin$Consequence=="synonymous_variant" & bin$sites=="Ti+Tv",]$avgRescaledHomRef)
        SYNdf$het= SYNdf$het +sum(bin[bin$Consequence=="synonymous_variant" & bin$sites=="Ti+Tv",]$avgRescaledHet)
        SYNdf$homAlt= SYNdf$homAlt +sum(bin[bin$Consequence=="synonymous_variant" & bin$sites=="Ti+Tv",]$avgRescaledHomAlt)
        # stop_gained:
        SGdf$homRef= SGdf$homRef + sum(bin[grepl("stop_gained",bin$Consequence) & bin$sites=="Ti+Tv",]$avgRescaledHomRef)
        SGdf$het= SGdf$het +sum(bin[grepl("stop_gained",bin$Consequence) & bin$sites=="Ti+Tv",]$avgRescaledHet)
        SGdf$homAlt= SGdf$homAlt +sum(bin[grepl("stop_gained",bin$Consequence) & bin$sites=="Ti+Tv",]$avgRescaledHomAlt)
        # transversions:
        MISdfTV$homRef= MISdfTV$homRef + sum(bin[grepl("missense_variant",bin$Consequence) & bin$sites=="TvOnly",]$avgRescaledHomRef)
        MISdfTV$het= MISdfTV$het +sum(bin[grepl("missense_variant",bin$Consequence) & bin$sites=="TvOnly",]$avgRescaledHet)
        MISdfTV$homAlt= MISdfTV$homAlt +sum(bin[grepl("missense_variant",bin$Consequence) & bin$sites=="TvOnly",]$avgRescaledHomAlt)
        # synonymous: ** note syn is treated differently to match with filter_vep -- no composite consequences allowed. but they are allowed for mis and sg*
        SYNdfTV$homRef= SYNdfTV$homRef + sum(bin[bin$Consequence=="synonymous_variant" & bin$sites=="TvOnly",]$avgRescaledHomRef)
        SYNdfTV$het= SYNdfTV$het +sum(bin[bin$Consequence=="synonymous_variant" & bin$sites=="TvOnly",]$avgRescaledHet)
        SYNdfTV$homAlt= SYNdfTV$homAlt +sum(bin[bin$Consequence=="synonymous_variant" & bin$sites=="TvOnly",]$avgRescaledHomAlt)
        # stop_gained:
        SGdfTV$homRef= SGdfTV$homRef + sum(bin[grepl("stop_gained",bin$Consequence) & bin$sites=="TvOnly",]$avgRescaledHomRef)
        SGdfTV$het= SGdfTV$het +sum(bin[grepl("stop_gained",bin$Consequence) & bin$sites=="TvOnly",]$avgRescaledHet)
        SGdfTV$homAlt= SGdfTV$homAlt +sum(bin[grepl("stop_gained",bin$Consequence) & bin$sites=="TvOnly",]$avgRescaledHomAlt)
        
      }}
      # update with final running totals
      MISdf$totalCDSSitesInBoot <- runningTotal
      SYNdf$totalCDSSitesInBoot <- runningTotal
      SGdf$totalCDSSitesInBoot <- runningTotal
      MISdfTV$totalCDSSitesInBoot <- runningTotal
      SYNdfTV$totalCDSSitesInBoot <- runningTotal
      SGdfTV$totalCDSSitesInBoot <- runningTotal
      
      GroupBoots=rbind(GroupBoots,MISdf,MISdfTV,SYNdf,SYNdfTV,SGdf,SGdfTV) # update GroupBoots
      
    }
    #GroupBoots$group <- group
    allBoots <- rbind(allBoots,GroupBoots) # add to allBoots
    write.table(GroupBoots,paste(out.dir,"/",outPREFIX,".",group,".allBoots.txt",sep=""),col.names = T,row.names = F,quote=F)
  }
write.table(allBoots,paste(out.dir,"/",outPREFIX,".Modern.Plus.Ancient.allBoots.txt",sep=""),col.names = T,row.names = F,quote=F)
