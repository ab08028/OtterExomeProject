require(ggplot2)
require(dplyr)
#require(grid)
#require(purrr) # for map, reduce
#require(readr) # for read_csv
#require(tidyr)# for unnest 
require(reshape2)
pops=c("CA","AK","AL","COM","KUR")
todaysdate=format(Sys.Date(),format="%Y%m%d")
genotypeDate=20181119
projection.date=20181221
hetFilter=0.5 # level that hets were filtered
plot.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/SFS/",genotypeDate,"/easySFS_projection/projection-",projection.date,"/hetFilter-",hetFilter,"/",sep="")
dir.create(plot.dir,recursive = T)
syn.data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/cds/synonymous/projection-",projection.date,"-hetFilter-",hetFilter,"/fastsimcoal2/",sep="")
mis.data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/cds/missense/projection-",projection.date,"-hetFilter-",hetFilter,"/fastsimcoal2/",sep="")
neut.data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/neutral/projection-",projection.date,"-hetFilter-",hetFilter,"/fastsimcoal2-1D-plusMonomorphic/",sep="")

# how do we want to plot? fsc format?


###### need to update these paths and plot all SFSes together (and try multiple filters) #######
################ Plot intermediate filtered SFSes #################
# ####################### fold sfs function ###############
# the SFSes are already folded from easySFS; but keep this function if you need to fold SFSes again sometime # 
# foldSFS <- function(sfs){
#   foldedSFS <- data.frame()
#   ss=length(sfs$frequency) - 1 # this is the ss in chromosomes
#   foldedBin=ss/2  # half as many as ss ; do seq from 0 --> foldedBins to include monomorphics
#   # note you have to start seq at 0 to include monomorphic bin 
#   for(i in seq(0,foldedBin)){
#     # if it's not the center point (which doesn't get added together)
#     # see wakeley coalescent eq 1.2
#     if(i==ss-i){
#       foldedSFS <- rbind.data.frame(foldedSFS, cbind.data.frame(frequency=i,count=(sfs[sfs$frequency==i,]$count + sfs[sfs$frequency==(ss-i),]$count)/2)) # if it's add midpoint (foldedLen, including 0 monomorphic bin), add together and divide by two (equivalent of not adding together)
#     }
#     else{
#       foldedSFS <- rbind.data.frame(foldedSFS, cbind.data.frame(frequency=i,count=(sfs[sfs$frequency==i,]$count + sfs[sfs$frequency==(ss-i),]$count))) # if not at mid point, just add together like normal
#     }
#   }
#   return(foldedSFS)
# }



########################### get all files ###############
# list of syn/mis files:
#syn.fileList=list.files(pattern=paste("pop0.obs",sep=""),path = syn.data.dir,full.names = F )
#mis.fileList=list.files(pattern=paste("pop0.obs",sep=""),path = mis.data.dir,full.names = F )
#neut.fileList=list.files(pattern=paste("pop0.obs",sep=""),path = neut.data.dir,full.names = F )

#"One limitation of the previous approach is that we don’t keep any auxilliary information #we may want to, such as the filenames of the files read. To keep the filename alongside #the data, we can read the data into a nested dataframe rather than a list, using the #mutate() function from dplyr. This gives us the following result:"
# https://serialmentor.com/blog/2016/6/13/reading-and-combining-many-tidy-data-files-in-R
# this assumes format is fsc folded format from EasySFS; if folded, there may be excess bins with 0s in them. dont' plot those; also don't plot monomorphic sites
for(pop in pops){
  
  ############## prepare the data ##############
  neut.input <- list.files(neut.data.dir,pattern=paste(pop,"_MAFpop0.obs",sep=""),full.names = T )
  syn.input <- list.files(syn.data.dir,pattern=paste(pop,"_MAFpop0.obs",sep=""),full.names = T )
  mis.input <- list.files(mis.data.dir,pattern=paste(pop,"_MAFpop0.obs",sep=""),full.names = T )
  neut.sfs <- melt(read.table(neut.input,skip = 1,header = T)) # skip first line: "13 folded "KUR""
  neut.sfs$label <- "neutral"
  syn.sfs <- melt(read.table(syn.input,skip = 1,header = T))
  syn.sfs$label <- "synonymous"
  mis.sfs <- melt(read.table(mis.input,skip = 1,header = T))
  mis.sfs$label <- "missense"
  # want to plot them separately and all together for the population
  # want to exclude d0_0 bin and 0-bins
  ########### put the syn, neut and mis sfses all together: ######
  
  all.sfs <- rbind.data.frame(syn.sfs,mis.sfs,neut.sfs)
  all.sfs$frequency <- sapply(strsplit(as.character(all.sfs$variable),"_"),"[",2)
  ######## exclude bins that are above the fold #########
  all.sfs_noMono_noZero <- all.sfs[all.sfs$value!=0 & all.sfs$frequency!="0",]
  ######## get proportional SFS as well: ########
  all.sfs_noMono_noZero_prop <- all.sfs_noMono_noZero %>%
    group_by(label) %>%
    mutate(proportion=value/sum(value))
  ####### order factors: #########
  all.sfs_noMono_noZero_prop$label <- factor(all.sfs_noMono_noZero_prop$label, levels=c("neutral","synonymous","missense"))
  ########### plot count-sfses all together ########
  sfsPlot <- ggplot(all.sfs_noMono_noZero_prop,aes(x=as.numeric(frequency),y=proportion,fill=label))+
    geom_bar(stat="identity",position = "dodge")+
    theme_bw()+
    theme(legend.title=element_blank(),legend.background = element_rect(fill = "transparent"),legend.position = c(0.5,0.9),legend.direction = "horizontal",legend.spacing.x = unit(.15,"cm"),text=element_text(size=18))+
    xlab("frequency")+
    ggtitle(paste(pop," folded SFS (projected)\nhet filter: ",hetFilter,sep=""))+
    scale_x_continuous(breaks=c(seq(1,max(as.numeric(all.sfs_noMono_noZero_prop$frequency)))))+
    scale_fill_manual(values=c("gray","dodgerblue","darkred"))
  sfsPlot
  ggsave(paste(plot.dir,pop,".cds.neutral.projected.prop.SFS.hetFilter.",hetFilter,".pdf",sep=""),sfsPlot,device = "pdf",height=5,width=7)
  
  ########## Plot each individually ############
  sfsPlot2 <- ggplot(all.sfs_noMono_noZero_prop,aes(x=as.numeric(frequency),y=value,fill=label))+
    geom_bar(stat="identity",position = "dodge")+
    theme_bw()+
    theme(legend.title=element_blank(),legend.background = element_rect(fill = "transparent"),legend.position = "none",legend.direction = "horizontal",legend.spacing.x = unit(.15,"cm"),text=element_text(size=18))+
    xlab("frequency")+
    ylab("count")+
    ggtitle(paste(pop," folded SFS (projected)\nhet filter: ",hetFilter,sep=""))+
    scale_x_continuous(breaks=c(seq(1,max(as.numeric(all.sfs_noMono_noZero_prop$frequency)))))+
    scale_fill_manual(values=c("gray","dodgerblue","darkred"))+
    facet_wrap(~label)
  sfsPlot2
  ggsave(paste(plot.dir,pop,".cds.neutral.projected.counts.faceted.SFS.hetFilter.",hetFilter,".pdf",sep=""),sfsPlot2,device = "pdf",height=5,width=12)
  
}

