############# 
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
              help="Number of bootstraps to perform per individual", metavar="numeric")
); 

# removing this option -- instead calculating driectly from point estimates based on bins:
#make_option(c("--avgSitesToDraw"), type="numeric", default=NULL, 
#            help="Number of sites to draw per individual (will be approximately this many, #not exactly) due to varying bin size. Based on average 'callable' cds sites calculated #elsewhere.", metavar="numeric")
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# test:
infile=opt$infile # output by step a-i - paste(out.dir,"/",outPREFIX,"Ind.",ind,".sumsPerBin.txt",sep="")
out.dir=opt$outdir
outPREFIX=opt$outPREFIX
numBoots=as.numeric(opt$numBoots) # do boots a number of times
#SitesToDraw=as.numeric(opt$avgSitesToDraw) getting this from point est now.

######### for testing #####
#allBins <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/bootStrapRegions/alt3test/test.Modern.Ancient.SumsPerIndividuals.RescaledPerWindow.binSize.1e+05.minIndPerWindow.9.minCallSitesPerWindow.500.txt",header=T)
#head(allBins)
#numBoots=5
#out.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/bootStrapRegions/alt3test"
#outPREFIX="test"                     

# read in results of step 7b-i:
allBins = read.table(infile,header=T,sep="\t")
# just pull out consequences we are interested in :
allBins_ConcOfInterest <- allBins[grepl("missense_variant",allBins$Consequence) | grepl("stop_gained",allBins$Consequence) | allBins$Consequence=="synonymous_variant",] # note missense and sg can be part of multiple annotations based on filter_vep (why I use grepl) but synonymous should just be alone (filter_vep), which is why we use == for it
allBins_ConcOfInterest$Consequence_BroadName <- as.character(allBins_ConcOfInterest$Consequence)
### rename composite consquences eg missense,splice_site just to be missense:
allBins_ConcOfInterest[grepl("missense_variant",allBins_ConcOfInterest$Consequence),]$Consequence_BroadName <- "missense_variant"
allBins_ConcOfInterest[grepl("stop_gained",allBins_ConcOfInterest$Consequence),]$Consequence_BroadName <- "stop_gained"
allBins_ConcOfInterest[allBins_ConcOfInterest$Consequence=="synonymous_variant",]$Consequence_BroadName <- "synonymous_variant"


###### get and write out point estimates from the bins: ###########
# okay get point estimates *per individual* (then average them to get modern/ancient average?)
pointEstimates <- allBins_ConcOfInterest %>% 
  group_by(group,sites,Consequence_BroadName,ind) %>%
  summarise(homRef=sum(sumHomRef_Rescaled),het=sum(sumHet_Rescaled),homAlt=sum(sumHomAlt_Rescaled)) # summing up over each category type that corresponds to missense,sg, syn. 
# get derived alleles:
pointEstimates$derivedAlleles <- (2*pointEstimates$homAlt) + pointEstimates$het
write.table(pointEstimates,paste(out.dir,"/",outPREFIX,".Modern.Ancient.PointEstimatesBasedonGoodBins.txt",sep=""),col.names = T,row.names = F,quote=F)
# also get callable sites in the bins passing filters (will be same across ancient and modern because you standardised bins)
callableSitesPerBin <- allBins_ConcOfInterest %>% 
  group_by(binNum) %>%
  summarise(unique(averageCalledSitesPerBin))
colnames(callableSitesPerBin) <- c("binNum","totalCallableSites")
totalCallableSitesAll <- sum(callableSitesPerBin$totalCallableSites)
SitesToDraw=totalCallableSitesAll # this is how many sites to draw for boots to match point ests
print(paste("Drawing ",SitesToDraw," sites for bootstraps",sep=""),quote=F)
##### randomly sample bins from an individual until some number is reached
allBootsPerInd=data.frame()
allBootsAvgs=data.frame()
print("starting resampling bootstraps",quote=F)
# do ancient and modern distributions separately
# want to pick the *same bins* for anc and modern -- not separately!!!

binNums=unique(allBins_ConcOfInterest$binNum)
for(i in seq(1,numBoots)){
  print(paste("starting bootstrap: ",i,sep=""),quote=F)
  runningTotal=0
  totalSitesInBoot=0
  numberOfDraws=0
  binsOfBoot=data.frame()
  ### start drawing bins:
  # this will slightly overshoot value, but that's okay, can rescale later if desired.
  while(runningTotal < SitesToDraw){
    # draw a random bin:
    randomBin=sample(binNums,1)
    numberOfDraws=numberOfDraws+1
    bin=allBins_ConcOfInterest[allBins_ConcOfInterest$binNum==randomBin,]
    # make sure what you're adding isn't going far over your SitesToDraw:
    # if adding the new amount won't make runningTotal > Sites to draw, then you can include that bin and update the running total
    runningTotal = runningTotal+unique(bin$averageCalledSitesPerBin) # this works; every line of df of bin has the total so you just need the unique value 
    # add the new bin to the rest of the bins:
    binsOfBoot=rbind(binsOfBoot,bin)
    # update running total:
    # get consequences of interest from the bin:
    totalSitesInBoot=runningTotal+totalSitesInBoot
  }
  # then once you have enough sites, summarise them all: 
  # test:
  # way more efficient: group by type of sites (Ti+Tv or TvOnly and by broad consequnce (mis, syn, sg)) and by ancient/modern
  # sum up avg GPs over each category: 
  sumsPerBoot <- binsOfBoot %>% 
    group_by(group,sites,Consequence_BroadName,ind) %>%
    summarise(homRef=sum(sumHomRef_Rescaled),het=sum(sumHet_Rescaled),homAlt=sum(sumHomAlt_Rescaled)) # summing up over each category type that corresponds to missense,sg, syn. 
  # need to get this per boot
  # get derived alleles:
  sumsPerBoot$derivedAlleles <- (2*sumsPerBoot$homAlt) + (sumsPerBoot$het)
  # update with final running totals
  sumsPerBoot$totalCDSSitesInBoot <- runningTotal
  # add metadata:
  sumsPerBoot$bootNum <- i
  sumsPerBoot$numberOfDraws <- numberOfDraws
  #head(sumsPerBoot)
  allBootsPerInd=rbind(allBootsPerInd,data.frame(sumsPerBoot)) # update GroupBoots
  ### Get averages for each set of three individuals for this particular bootstrap
  avgsPerGroupPerBoot <- sumsPerBoot %>%
    group_by(group,sites,Consequence_BroadName,totalCDSSitesInBoot) %>%
    summarise(avgHomRef=mean(homRef),avgHet=mean(het),avgHomAlt=mean(homAlt),avgDerivedAlleles=mean(derivedAlleles))
  avgsPerGroupPerBoot$bootNum <- i
  avgsPerGroupPerBoot$numberOfDraws <- numberOfDraws
  allBootsAvgs = rbind(allBootsAvgs,data.frame(avgsPerGroupPerBoot))
}

## add last bits of metadata:
allBootsPerInd$SitesToDraw <- SitesToDraw
allBootsPerInd$derivedAlleles_Rescaled <- (allBootsPerInd$derivedAlleles/allBootsPerInd$totalCDSSitesInBoot) * allBootsPerInd$SitesToDraw # rescale by SitesToDraw to get spot on estimates of exaclty the SitesToDRaw amount per bootstrap
allBootsPerInd$homAlt_Rescaled <- (allBootsPerInd$homAlt/allBootsPerInd$totalCDSSitesInBoot) * allBootsPerInd$SitesToDraw  # rescale by exact amnt of desired bootstrap sites.

allBootsAvgs$SitesToDraw <- SitesToDraw
# rescale by exact estimate of how many sites drawn in the bootstrap
allBootsAvgs$derivedAlleles_Rescaled <- (allBootsAvgs$avgDerivedAlleles/allBootsAvgs$totalCDSSitesInBoot) * allBootsAvgs$SitesToDraw
# homAlt rescaled by 
allBootsAvgs$homAlt_Rescaled <- (allBootsAvgs$avgHomAlt/allBootsAvgs$totalCDSSitesInBoot) * allBootsAvgs$SitesToDraw
  
# write out:
write.table(allBootsPerInd,paste(out.dir,"/",outPREFIX,".perIndividualBootstrapSums.txt",sep=""),col.names = T,row.names = F,quote=F)
# write out:
write.table(allBootsAvgs,paste(out.dir,"/",outPREFIX,".perGroupAveragesBootstraps.useThis.txt",sep=""),col.names = T,row.names = F,quote=F)
