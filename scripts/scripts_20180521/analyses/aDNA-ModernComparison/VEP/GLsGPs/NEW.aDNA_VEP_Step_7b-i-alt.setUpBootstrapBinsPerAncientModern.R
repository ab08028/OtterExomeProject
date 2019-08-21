#BiocManager::install("plyranges")
require(GenomicRanges)
require(dplyr)
require(plyranges)
library("optparse")
option_list = list(
  make_option(c("--infile"), type="character", default=NULL, 
              help="path to your input file (superfile of cds sites with annotations from vep and GPs from angsd as well as depths", metavar="file"),
  make_option(c("--chrSizes"), type="character", default=NULL, 
              help="path to file with mustela chromosome sizes", metavar="character"),
  make_option(c("--outdir"), type="character", default=NULL, 
              help="path to outdir", metavar="dir"),
  make_option(c("--outPREFIX"), type="character", default=NULL, 
              help="outfilePrefix", metavar="prefix"),
  make_option(c("--minDepth"), type="numeric", default=NULL, 
              help="min depth per individual for a site to be 'callalble'", metavar="numeric"),
  make_option(c("--minGP"), type="numeric", default=NULL, 
              help="minimum value of the max. GP per site, per individual. Use 0.95.", metavar="numeric"),
  make_option(c("--binsize"), type="numeric", default=NULL, 
              help="Size of bin to chunk the genome into (should be > than a recombination block)", metavar="numeric"),
  #make_option(c("--indNum"), type="numeric", default=NULL, 
  #             help="Individual number assigned by ANGSD, starts at 0 (for my study it's 0-8))", metavar="numeric"),
  make_option(c("--bamList"), type="character", default=NULL, 
              help="path to list of bams in angsd (gives IDs) ***ASSUMES THAT ANCIENT SAMPLE IDs start with 'A'!!!!!!#", metavar="file")
); 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cdsFile=opt$infile
mustelaChrSizesFile=opt$chrSizes
out.dir=opt$outdir
outPREFIX=opt$outPREFIX
minDepth=as.numeric(opt$minDepth) # make this match whatever I used to get point estimate
minGP=as.numeric(opt$minGP) # make this match whatever I used to get point estimate
binsize=as.numeric(opt$binsize) # start with 100kb
ind=as.numeric(opt$indNum)
bamListFile=as.character(opt$bamList)

transversions=c('A,C','C,A','A,T','T,A','C,G','G,C','G,T','T,G')

# for testing:
#minDepth=2 # make this match whatever I used to get point estimate
#minGP=0.95 # make this match whatever I used to get point estimate
#binsize=100000 # start with 100kb
#bamListFile="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/SampleIDsInOrder.LowCoverageOnly.BeCarefulOfOrder.txt"
bamList=read.table(bamListFile)
# need to turn bamList into something with indexes and get the modern and ancient indices from it. 
colnames(bamList) <- "sample"
# add numbers (0 based)
bamList$indNumber <- unlist(seq(0,length(bamList$sample)-1)) # from angsd
# add group info **** THIS ASSUMES THAT ANCIENT SAMPLES START WITH ^A***
bamList$group <- "Modern"
bamList[grepl("^A",bamList$sample),]$group <- "Ancient"
# get lists of IDs that correspond to anc and mod:
ancientIDs=bamList[bamList$group=="Ancient",]$indNumber
modernIDs=bamList[bamList$group=="Modern",]$indNumber

print("WARNING: NOTE THAT THIS SCRIPT ASSUMES aDNA IDs start with 'A' -- this is specific to my dataset and my not be specific to others!!!!!!" )
print(paste("ancient samples are flagged as:", bamList[bamList$group=="Ancient",]$sample))
print(paste("modern samples are flagged as: ",bamList[bamList$group=="Modern",]$sample))
print("If this isn't right, you must modify script -- sorry in advance")
###############read in cds superfile ######################
#cds=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/sandbox/bootStrapRegions/test.out.bed",header=T,sep="\t",strip.white = T, comment.char = "") # converts # to X. in table ### will R be able to read it in on Hoffman?
cds=read.table(cdsFile,header=T,sep="\t",strip.white = T, comment.char = "")
# pull out the variables we care out (all Ind information, site information and Consequences from VEP):
Variables=c("chromo","start0based","end","position","major","minor","ref",names(cds)[grep("^Ind",names(cds))],names(cds)[grep("TotDepth",names(cds))],"Consequence","Extra") # note: IndX is 00, IndX.1 is 01, IndX.2 is 11 
# just pull those columsn so cds is smaller:
cds2=cds[,Variables]
## make cds file a Granges object:
cds_GRanges <- makeGRangesFromDataFrame(df=cds2,keep.extra.columns = T,ignore.strand = T,seqnames.field = "chromo",start.field = "start0based",end.field = "end",starts.in.df.are.0based = T,seqinfo = NULL)
############# get cds superfile and make it a Granges object #######
mustelaChrSizes=read.table(mustelaChrSizesFile,sep=",",stringsAsFactors = F,strip.white = T)
#mustelaChrSizes <- (read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/SlidingWindowheterozygosity/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.224.ChrLengths.txt",sep=",",stringsAsFactors = F,strip.white = T))
mustelaChrSizes <- data.frame(t(mustelaChrSizes))
colnames(mustelaChrSizes) <- "record"
mustelaChrSizes$scaff <- lapply(strsplit(as.character(mustelaChrSizes$record),":"),"[",1)
mustelaChrSizes$size <- lapply(strsplit(as.character(mustelaChrSizes$record),":"),"[",2)

####### intersect with the scaffolds that are in my cds Granges object #########
mustelaChrSizes_cds <-  mustelaChrSizes[mustelaChrSizes$scaff %in% unique(cds$chromo),] # 
sizes <- as.numeric(mustelaChrSizes_cds$size)
names(sizes) <- mustelaChrSizes_cds$scaff
########### split scaffolds into even sized bins (some will be empty because don't contain cds) #######
###  bins for bootstrapping
bins   <- tileGenome(sizes, tilewidth=binsize, cut.last.tile.in.chrom=T)
#bins
# add a bin number to each bin
bins$binNum <- seq(1,length(bins))
print(c("number of unique bins:", length(unique(bins$binNum))))
##### UPDATE:
write.table(bins,paste(out.dir,"/",outPREFIX,".BinCoords.binSize.",binsize,".txt",sep=""),col.names = T,row.names = F,quote=F,sep="\t")
############### want to sum stuff up per bin per category #######
##### loop over bins (make run in parallel?) ####
allBinsallInds=data.frame()
for(bin in seq(1,length(unique(bins$binNum)))){
#for(bin in seq(1,10)){
  print(paste("starting bin ",bin),quote=F)
  binTotals=data.frame()  # reset this for every bin
  callableSiteTotalsPerBin=data.frame() # just for totals per individual
  # go through each individual
  subset <- subsetByOverlaps(cds_GRanges, bins[bins$binNum==bin,])
  # go through inds:
  for(ind in seq(0,length(bamList$sample)-1)){
    #print(paste("starting ind ",ind),quote=F)
    indTotals=data.frame() # per ind
    print(paste("starting ind",ind))
    IndVariables=c("seqnames","start","end","position","major","minor","ref",paste("Ind",ind,sep=""),paste("Ind",ind,".1",sep=""),paste("Ind",ind,".2",sep=""),paste("ind",ind,"TotDepth",sep=""),"Consequence","Extra") # note: IndX is 00, IndX.1 is 01, IndX.2 is 11 
    indOnly <- data.frame(subset)[,IndVariables] # . 
    # make generic column names:
    colnames(indOnly) <- c("chromo","start0based","end","position1based","major","minor","ref","homRef","het","homAlt","indDepth","Consequence","Extra")
    ####### employ filters:
    # get maxGP value among the homRef, het or homAlt GTs:
    indOnly$maxGP <- apply(indOnly[,c("homRef","het","homAlt")],1,max) # this works
    # 1. min depth
    # 20180816: fixed this to be >= not just > -- this could cause discrepancy! 
    # also checks for canonical (just added this --20190820)
    indOnly_filter <- indOnly %>% 
      filter(indDepth >= minDepth) %>%
      filter(maxGP >= minGP) %>% 
      filter(grepl("CANONICAL=YES",Extra))
    # this makes sure that whatever the largest GP is for the site, that it is > some cutoff for the site (makes a 'callable site')
    # add canonical filter
    
    # okay so now it's been filtered by GP and by Depth
    # want to sum up by category
    totalCallableSites=dim(indOnly_filter)[1] # after filtering, the sites remaining are the total callable sites. 
    # Exclude individual if they don't have any callable sites in the window (saves NAs)
    if(totalCallableSites>0){
      # add to callable sites df:
      callableSiteTotalsPerBin = rbind(callableSiteTotalsPerBin,data.frame(total=totalCallableSites,bin=bin,ind=ind))
      # make sure ref is in major minor
      indOnly_filter$Alleles <- paste(indOnly_filter$major,indOnly_filter$minor,sep=",")
      #if(totalCallableSites>0){
      # all sites
      # use .drop=F in group_by to get zeros if it's an empty category -- helps avoid empty dfs if sites are missing; will put in zero counts instead
      indTotals <- indOnly_filter %>% 
        group_by(Consequence,.drop=F) %>%
        summarise(sumHomRef=sum(homRef),sumHet=sum(het),sumHomAlt=sum(homAlt))
      indTotals$ind <- ind
      indTotals$sites <- "Ti+Tv"
      indTotals$totalCallableSitesPerBin <- totalCallableSites
      binTotals <- rbind(binTotals,indTotals)
      
      # transverions:
      indTotalsTV <- indOnly_filter %>% 
        group_by(Consequence,.drop=F) %>%
        filter(Alleles %in% transversions) %>%
        summarise(sumHomRef=sum(homRef),sumHet=sum(het),sumHomAlt=sum(homAlt))
      indTotalsTV$ind <- ind
      indTotalsTV$sites <- "TvOnly"
      indTotalsTV$totalCallableSitesPerBin <- totalCallableSites
      binTotals <- rbind(binTotals,indTotalsTV)
      
    }}
  # want to count up individuals 
  ##### need to do some normalizing ####
  # get avg sites per group of individuals:
  # count up how many inds have data (so can filter on later if you want)
  totalIndsWithData=length(unique(binTotals$ind)) # counts up inds with data. If no individuals have data, skip the bin:
  if(totalIndsWithData>0){
    # add that info to df
    # divide each by total callable sites per individual
    # generates NaN if totalCallableSitesPerBin is 0
    binTotals$sumHomRef_Frac <- binTotals$sumHomRef/binTotals$totalCallableSitesPerBin
    binTotals$sumHet_Frac <- binTotals$sumHet/binTotals$totalCallableSitesPerBin    
    binTotals$sumHomAlt_Frac <- binTotals$sumHomAlt/binTotals$totalCallableSitesPerBin
    # then multiply by average sites in the bin for modern or for ancient (for now keeping separate)
    ### get averages per bin: (useing separate callableSiteTotalsPerBin df so that multiple entries don't get counted, just one total per individual for mean )
    ########## get modern and ancient average called sites per bin (doing separately for now) ####
    # update for modern and ancient 
    averageModernPerBin=mean(callableSiteTotalsPerBin[callableSiteTotalsPerBin$ind %in% modernIDs,]$total,na.rm=T)
    # ancient:
    averageAncientPerBin=mean(callableSiteTotalsPerBin[callableSiteTotalsPerBin$ind %in% ancientIDs,]$total,na.rm=T)
    # add avg sites info to binTotals
    binTotals$averageCalledSitesPerBin <- NA
    binTotals[binTotals$ind %in% modernIDs,]$averageCalledSitesPerBin <- averageModernPerBin
    # add avg sites info to binTotals
    binTotals[binTotals$ind %in% ancientIDs,]$averageCalledSitesPerBin <- averageAncientPerBin
    
    # then want to rescale values:
    binTotals$sumHomRef_Rescaled <- binTotals$sumHomRef_Frac * binTotals$averageCalledSitesPerBin
    binTotals$sumHet_Rescaled <- binTotals$sumHet_Frac * binTotals$averageCalledSitesPerBin
    binTotals$sumHomAlt_Rescaled <- binTotals$sumHomAlt_Frac * binTotals$averageCalledSitesPerBin
    
    ### add individual group info ###
    binTotals$group <- NA
    binTotals[binTotals$ind %in% modernIDs,]$group <- "Modern"
    binTotals[binTotals$ind %in% ancientIDs,]$group <- "Ancient"
    ## need to then average over those:
    summaries <- binTotals %>% 
      group_by(group,sites,Consequence,averageCalledSitesPerBin,.drop=F) %>%
      summarise(avgRescaledHomRef=mean(sumHomRef_Rescaled,na.rm=T),avgRescaledHet=mean(sumHet_Rescaled,na.rm=T),avgRescaledHomAlt=mean(sumHomAlt_Rescaled,na.rm=T))
    # checked this, it works. cool.
    ### add metadata:
    summaries$binNum <- bin
    summaries$scaff <- as.character(unlist(seqnames(bins[bins$binNum==bin,])))
    summaries$start <- start(bins[bins$binNum==bin,])
    summaries$end <- end(bins[bins$binNum==bin,])
    summaries$width <- width(bins[bins$binNum==bin,])
    summaries$minDepth <- minDepth
    summaries$minGP <- minGP
    summaries$totalIndsWithData <- totalIndsWithData
    # finally: add to allBinsallInds df:
    allBinsallInds <- rbind(allBinsallInds,data.frame(summaries))
    
    
  }}
# write out modern and ancient:
# could do NA.omit? 
write.table(allBinsallInds,paste(out.dir,"/",outPREFIX,".Modern.Ancient.AvgsPerGroup.PerBinbinSize.",binsize,".txt",sep=""),col.names = T,row.names = F,quote=F,sep="\t")

write.table(allBinsallInds[allBinsallInds$group=="Modern",],paste(out.dir,"/",outPREFIX,".ModernOnly.AvgsPerGroup.PerBin.binSize.",binsize,".txt",sep=""),col.names = T,row.names = F,quote=F,sep="\t")

write.table(allBinsallInds[allBinsallInds$group=="Ancient",],paste(out.dir,"/",outPREFIX,".AncientOnly.AvgsPerGroup.PerBin.binSize.",binsize,".txt",sep=""),col.names = T,row.names = F,quote=F,sep="\t")
# what needs to be written out?

