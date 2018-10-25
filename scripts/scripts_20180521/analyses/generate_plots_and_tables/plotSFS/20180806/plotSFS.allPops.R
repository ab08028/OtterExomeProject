### Plot the SFS
require(ggplot2)
require(dplyr)
require(gridExtra)
# 20181019 updating this script to work from output of *my* SFS script, not the Tanya script. These SFSs were made from these neutral VCF files : ${pop}_neutral_all_9_rmAllHet_rmRelativesAdmixed_passingAllFilters_allCalled.vcf.gz

pops=c("CA","AK","AL","COM","KUR")  # your populations
rundate=20180806 # genotype calling date
sfsdate=20181019 # date sfses were made
suffix=paste("sfs.R.format.",sfsdate,".txt",sep="")
############### set up your colors -- keep this consistent across all plots ######
require(RColorBrewer)
colorPal=RColorBrewer::brewer.pal(n=length(pops),name = "Dark2")
colors=list(CA=colorPal[1],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors
#######################################################################################
# set up dir structure 
todaysdate=format(Sys.Date(),format="%Y%m%d") # date you make plots

# post filtering and removing admixed invididuals etc.:
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",rundate,"/neutralSFS/",sep="")

plot.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/SFS/",rundate,"/neutralSFS/",sep="")

dir.create(plot.dir,recursive = T,showWarnings = F)
dir.create(int.plot.dir,recursive = T,showWarnings = F)


# make a little function to fold SFS that is set up as frequency count (I checked and it matches dadi folded sfs)
foldSFS <- function(sfs){
  foldedSFS <- data.frame()
  ss=length(sfs$frequency) - 1 # this is the ss in chromosomes
  foldedBin=ss/2  # half as many as ss ; do seq from 0 --> foldedBins to include monomorphics
  # note you have to start seq at 0 to include monomorphic bin 
  for(i in seq(0,foldedBin)){
    # if it's not the center point (which doesn't get added together)
    # see wakeley coalescent eq 1.2
    if(i==ss-i){
      foldedSFS <- rbind.data.frame(foldedSFS, cbind.data.frame(frequency=i,count=(sfs[sfs$frequency==i,]$count + sfs[sfs$frequency==(ss-i),]$count)/2)) # if it's add midpoint (foldedLen, including 0 monomorphic bin), add together and divide by two (equivalent of not adding together)
    }
    else{
      foldedSFS <- rbind.data.frame(foldedSFS, cbind.data.frame(frequency=i,count=(sfs[sfs$frequency==i,]$count + sfs[sfs$frequency==(ss-i),]$count))) # if not ad mid point, just add together like normal
    }
    }
  return(foldedSFS)
}


#### Plot all SFSes: 
allSFS_list=list()
allPlots_list=list()


# this is a nice for-loop; it generates and writes out a plot for each population
for(i in (1:length(pops))){
  pop=pops[i]
  sfs <- read.table(paste(data.dir,pop,".unfolded.",suffix,sep=""),header=T,stringsAsFactors = F)
  sfs_folded <- foldSFS(sfs)
  # write out folded SFS:
  ss=dim(sfs)[1]-1 # sample size in chromsomes (haps). full size of sfs is ss+1 because of 0 bin
  # exclude monomorphic at freq 0 and total bins -1: 
  sfsNoMono <- sfs_folded[sfs_folded$frequency!=0 & sfs_folded$frequency!=ss,]
  sfsNoMono$population <- pop
  popcolor=colors[pop]
  allSFS_list[[i]] <- sfsNoMono # add the SFS to your list
  p <- ggplot(allSFS_list[[i]],aes(x=as.numeric(frequency),y=count))+
    geom_bar(stat="identity",fill=popcolor)+
    scale_x_continuous(breaks=c(seq(1,nrow(allSFS_list[[i]]))))+
    ggtitle(paste(pop,": Folded Neutral SFS",sep=""))+
    theme_bw()+
    facet_wrap(~population,scales="free_x")+
    xlab("frequency")+
    ylab("SNP Count")
  allPlots_list[[i]]=p
  ggsave(paste(plot.dir,pop,"_foldedNeutralSFS_counts.",todaysdate,".pdf",sep=""),p,device="pdf",width=7,height=5)
}
allPlots_list
# can put all together into a dataframe:
allSFS <- bind_rows(allSFS_list)

################################# Proportional SfS ###################

allSFS_list=list()
allPlots_list=list()


# this is a nice for-loop; it generates and writes out a plot for each population
for(i in (1:length(pops))){
  pop=pops[i]
  sfs <- read.table(paste(data.dir,pop,".unfolded.",suffix,sep=""),header=T,stringsAsFactors = F)
  sfs_folded <- foldSFS(sfs)
  ss=dim(sfs)[1]-1 # sample size in chromsomes (haps). full size of sfs is ss+1 because of 0 bin
  # exclude monomorphic at freq 0 and total bins -1: 
  sfsNoMono <- sfs_folded[sfs_folded$frequency!=0 & sfs_folded$frequency!=ss,]
  sfsNoMono$population <- pop
  sfsNoMono$proportion <- sfsNoMono$count/sum(sfsNoMono$count)
  write.table(sfsNoMono[,c("frequency","proportion")],paste(data.dir,pop,".folded.PROPORTIONS.",suffix,sep=""),quote=F,row.names = F,sep="\t")
  popcolor=colors[pop]
  allSFS_list[[i]] <- sfsNoMono # add the SFS to your list
  p <- ggplot(allSFS_list[[i]],aes(x=as.numeric(frequency),y=proportion))+
    geom_bar(stat="identity",fill=popcolor)+
    scale_x_continuous(breaks=c(seq(1,nrow(allSFS_list[[i]]))))+
    ggtitle(paste(pop,": Folded Neutral SFS",sep=""))+
    theme_bw()+
    facet_wrap(~population,scales="free_x")+
    xlab("frequency")+
    ylab("Proportion")
  allPlots_list[[i]]=p
  ggsave(paste(plot.dir,pop,"_foldedNeutralSFS_proportions.",todaysdate,".pdf",sep=""),p,device="pdf",width=7,height=5)
}
allPlots_list
# can put all together into a dataframe:
allSFS <- bind_rows(allSFS_list)


# get total SNPs:
snpCount <- allSFS %>% 
  group_by(population) %>%
  tally(sum(count)) 
sSize <- allSFS %>% 
  group_by(population) %>%
  tally() 
########## admixed SFS ############
# ADprefix="all_9_rmAllHet_passingAllFilters_allCalled"
# ADpops=c("AK","KUR") # admixed pops
# ADallSFS_list=list()
# ADallPlots_list=list()
# for(i in (1:length(ADpops))){
#   pop=ADpops[i]
#   sfs <- read.table(paste(AD.data.dir,"/admixIndOnly_",pop,"_",ADprefix,".filtered.sfs.",sfsdate,".out",sep=""),header=T,stringsAsFactors = F) 
#   sfs$population <- pop
#   popcolor=colors[pop]
#   ADallSFS_list[[i]] <- sfs # add the SFS to your list
#   p <- ggplot(ADallSFS_list[[i]],aes(x=as.numeric(frequency),y=num_variants))+
#     geom_bar(stat="identity",fill=popcolor)+
#     scale_x_continuous(breaks=c(seq(1,nrow(ADallSFS_list[[i]]))))+
#     ggtitle(paste("Admixed Individuals from ",pop,": Folded Neutral SFS",sep=""))+
#     theme_bw()+
#     facet_wrap(~population,scales="free_x")+
#     xlab("frequency")+
#     ylab("SNP Count")
#   ADallPlots_list[[i]]=p
#   ggsave(paste(plot.dir,"/admixed/",pop,"_foldedNeutralSFS_counts.admixedIndividuals.",todaysdate,".pdf",sep=""),p,device="pdf",width=7,height=5)
# }
########## sandbox: intermediate sfses ############
######## Look at intermediate steps:
# /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/SFS/20180806/neutralSFS/filteringSteps
# use this for intermediate sfs: sfs <- read.table(paste(int.data.dir,pop,"_all_7_passingAllFilters_allCalledraw_variants.sfs.out",sep=""),header=T,stringsAsFactors = F)
# allSFS_list=list()
# allPlots_list=list()
# ###### something is odd here -- this is not for neutral regions; is for whole capture region; for some reason the even frequencies are way inflated.
# # I think it may be because tanya's script wants all sites called, so it won't know what to do with no call genotypes?
# for(i in (1:length(pops))){
#   pop=pops[i]
#   sfs1 <- read.table(paste(int.data.dir,pop,"_all_1_TrimAlt_raw_variants.sfs.out",sep=""),header=T,stringsAsFactors = F) 
#   sfs1$population <- pop
#   sfs1$label <- "1. raw"
#   sfs5 <- read.table(paste(int.data.dir,pop,"_all_5_passingFilters_raw_variants.sfs.out",sep=""),header=T,stringsAsFactors = F) 
#   sfs5$population <- pop
#   sfs5$label <- "5. passing initial filters (old scheme)"
#   sfs6 <- read.table(paste(int.data.dir,pop,"_all_6_rmBadIndividuals_passingFilters_raw_variants.sfs.out",sep=""),header=T,stringsAsFactors = F)
#   sfs6$population <- pop
#   sfs6$label <- "6. remove low coverage individuals"
#   sfs7 <- read.table(paste(int.data.dir,pop,"_all_7_passingBespoke_maxNoCallFrac_0.2_rmBadIndividuals_passingFilters_raw_variants.sfs.out",sep=""),header=T,stringsAsFactors = F)
#   sfs7$population <- pop
#   sfs7$label <- "7. bespoke filters (still has admixed/relatives/all-het sites)"
#   popcolor=colors[pop]
#   # rbind all those sfses
#   sfs <- rbind(sfs1,sfs5,sfs6,sfs7)
#   allSFS_list[[i]] <- sfs # add the set of SFSes to your list
#   p <- ggplot(allSFS_list[[i]],aes(x=as.numeric(frequency),y=num_variants,fill=as.factor(label)))+
#     geom_bar(stat="identity",position="dodge")+
#     scale_x_continuous(breaks=c(seq(1,nrow(allSFS_list[[i]]))))+
#     ggtitle(paste(pop,": Folded Neutral SFS",sep=""))+
#     theme_bw()+
#     facet_wrap(~label,scales="free_x")+
#     xlab("frequency")+
#     ylab("SNP Count")
#   allPlots_list[[i]]=p
# }
#   #ggsave(paste(plot.dir,pop,"_foldedNeutralSFS_counts.",todaysdate,".pdf",sep=""),p,device="pdf",width=7,height=5)
# 
# # can also arrange as a grid:
# 
# # do.call("grid.arrange",c(allPlots_list))
# 
# # can also plot as facetted:
# #allSFS <- bind_rows(allSFS_list)
# 
# #p1 <- ggplot(allSFS,aes(x=as.numeric(frequency),y=num_variants,fill=population))+
#   #geom_bar(stat="identity")+
#  #scale_x_continuous(breaks=c(1,nrow(allSFS)))+
#   #ggtitle(paste("Folded Neutral SFS",sep=""))+
#   #theme_bw()+
#   #facet_wrap(~population,scales="free_x")+
#   #xlab("frequency")+
#   #ylab("SNP Count")
# #p1
