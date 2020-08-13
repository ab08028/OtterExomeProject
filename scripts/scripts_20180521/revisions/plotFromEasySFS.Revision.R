require(reshape2)
require(ggplot2)
require(RColorBrewer)
############# Plot **NEUTRAL **FSC format SFS in R ############
############### set up your colors -- keep this consistent across all plots ######
require(RColorBrewer)
colorPal=RColorBrewer::brewer.pal(n=8,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[7],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors

######### file specifics ######
datadirs=c("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/SFS_BasedOnSSONeutralRegions/20200719_SSO/easySFS/neutral/projection-20200729-hetFilter-0.75/","/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/")

# these are combos of gentoype dates (first) and projection dates (second)

#filter="hetFilter-0.75" # date you did the projection (so you can try different projections)
plot.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/SFS_BasedOnSSONeutralRegions/plots",sep="")
dir.create(plot.dir,recursive = T)
popFile=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/easySFSPopMapFiles/samplesPop.Headers.forEasySFS.3.20181119.txt",header=F) # used in easySFS; update if used a different one. The order of populations in this script is the order that pops are assigned numbers in easy sfs (that's a bit hacky). So get order of pops from this file for 0/1 (double check manually)
colnames(popFile) <- c("sample","population")
popOrder <- as.character(unique(popFile$population)) # this should be CA,AK,AL,COM,KUR  for my project
allSFSes <- data.frame()
for(data.dir in datadirs){
#data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/neutral/projection-",dateOfProjection,"-",filter,"/",sep="")
fsc.format.dir=paste(data.dir,"fastsimcoal2/",sep="")
dadi.format.dir = paste(data.dir,"dadi/",sep="")

###################### fold sfs function ############
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
      foldedSFS <- rbind.data.frame(foldedSFS, cbind.data.frame(frequency=i,count=(sfs[sfs$frequency==i,]$count + sfs[sfs$frequency==(ss-i),]$count))) # if not at mid point, just add together like normal
    }
  }
  return(foldedSFS)
}
############ Plot 1D sfses for each population ###############################
# melt sfs:
for(pop in popOrder) {
  input <- list.files(fsc.format.dir,pattern=pop,full.names = T)
  sfs <- read.table(input,skip = 1,header = T) # skip first line: "1 observation"
  sfs <- melt(sfs) # melt it to make it longways
  # get frequency from the d0_0 
  sfs$frequency <- as.numeric(sapply(strsplit(as.character(sfs$variable),"_"),"[",2))
  sfs <- sfs[,c("frequency","value")]
  colnames(sfs) <- c("frequency","count") # so it will work in the foldSFS function
  sfs_fold <- foldSFS(sfs) # gets rid of the trailing zeroes (was folded already but wanted to elimiante extra zeroes)
  sfs_exclMono <- sfs_fold[sfs_fold$frequency!=0,]
  sfs_exclMono$datadir <- data.dir
  sfs_exclMono$pop <- pop
  sfs_exclMono$proportion <- sfs_exclMono$count/sum(sfs_exclMono$count)
  # p <- ggplot(sfs_exclMono,aes(x=frequency, y=count,fill=pop))+
  #   geom_bar(stat="identity")+
  #   theme_bw()+
  #   scale_x_continuous(breaks=c(seq(1,length(sfs_exclMono$frequency))))+
  #   ggtitle(pop) +
  #   scale_fill_manual(values=colors[pop])
  allSFSes <- rbind(allSFSes,sfs_exclMono)
 # ggsave(paste(plot.dir,pop,"folded.easySFS.proj.sfs.pdf",sep=""),device="pdf",height=5,width=7)
}
}
########### NEED PROPORTIONAL SFS!!!
## add a label:
allSFSes$ref <- NA
allSFSes[allSFSes$datadir=="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/",]$ref <- "ferret" 
allSFSes[allSFSes$datadir=="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/revisions/SFS_BasedOnSSONeutralRegions/20200719_SSO/easySFS/neutral/projection-20200729-hetFilter-0.75/",]$ref <- "sea otter" 
for(pop in popOrder){
  p1 <- ggplot(allSFSes[allSFSes$pop==pop,],aes(x=frequency, y=proportion,fill=ref))+
    geom_bar(stat="identity",position="dodge")+
    theme_classic()+
    scale_x_continuous(breaks=c(seq(1,length(allSFSes[allSFSes$pop==pop,]$frequency))))+
    ggtitle(pop)+
    theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),text = element_text(size=20),axis.text =element_text(size=20))
  p1
  ggsave(paste(plot.dir,pop,".compare.Mfur.SSO.Refs.pdf",sep=""),p1,device="pdf",height=5,width=7)
  
  p2 <- ggplot(allSFSes[allSFSes$pop==pop,],aes(x=frequency, y=count,fill=ref))+
    geom_bar(stat="identity",position="dodge")+
    theme_classic()+
    scale_x_continuous(breaks=c(seq(1,length(allSFSes[allSFSes$pop==pop,]$frequency))))+
    ggtitle(pop)+
    theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),text = element_text(size=20),axis.text =element_text(size=20))
  p2
  ggsave(paste(plot.dir,pop,".compare.Mfur.SSO.Refs.COUNTS.pdf",sep=""),p2,device="pdf",height=5,width=7)
  
}

