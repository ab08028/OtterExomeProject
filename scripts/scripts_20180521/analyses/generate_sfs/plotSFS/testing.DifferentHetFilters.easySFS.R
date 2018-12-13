####### troubleshooting het filters #########

# 80% het filter:
require(reshape2)
require(ggplot2)
require(RColorBrewer)
############# Plot FSC format SFS in R ############
############### set up your colors -- keep this consistent across all plots ######
require(RColorBrewer)
colorPal=RColorBrewer::brewer.pal(n=6,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[7],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors

######### file specifics ######
genotypeDate="20181119"
dateOfProjection="20181212" # date you did the projection (so you can try different projections)
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/troubleshoot-hetFilter/",sep="")
plot.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/SFS/",genotypeDate,"/easySFS_projection/troubleshoot-hetFilter/",sep="")
dir.create(plot.dir,recursive = T)
fsc.format.dir=paste(data.dir,"fastsimcoal2/",sep="")
dadi.format.dir = paste(data.dir,"dadi/",sep="")
popFile=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/information/samples/easySFSPopMapFiles/samplesPop.Headers.forEasySFS.3.20181119.txt",header=F) # used in easySFS; update if used a different one. The order of populations in this script is the order that pops are assigned numbers in easy sfs (that's a bit hacky). So get order of pops from this file for 0/1 (double check manually)
colnames(popFile) <- c("sample","population")
popOrder <- as.character(unique(popFile$population)) # this should be CA,AK,AL,COM,KUR  for my project

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
  p <- ggplot(sfs_exclMono,aes(x=frequency, y=count,fill=pop))+
    geom_bar(stat="identity")+
    theme_bw()+
    scale_x_continuous(breaks=c(seq(1,length(sfs_exclMono$frequency))))+
    ggtitle(pop)+
    scale_fill_manual(values=colors[pop])
  ggsave(paste(plot.dir,pop,"folded.easySFS.proj.sfs.pdf",sep=""),device="pdf",height=5,width=7)
}


############### plot all 3 next to each other #########
sfs1 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/projection-20181212/fastsimcoal2/COM_MAFpop0.obs",skip=1,header=T)
sfs0.8 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/troubleshoot-hetFilter-0.8/fastsimcoal2/COM_MAFpop0.obs",skip = 1,header = T)
sfs0.9 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/troubleshoot-hetFilter-0.9/fastsimcoal2/COM_MAFpop0.obs",skip = 1,header = T)
sfs0.7 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/troubleshoot-hetFilter-0.7/fastsimcoal2/COM_MAFpop0.obs",skip = 1,header = T)
sfs1 <- melt(sfs1) # melt it to make it longways
sfs0.8 <- melt(sfs0.8)
sfs0.9 <- melt(sfs0.9)
sfs0.7 <- melt(sfs0.7)
# get frequency from the d0_0 
sfs1$frequency <- as.numeric(sapply(strsplit(as.character(sfs1$variable),"_"),"[",2))
sfs1 <- sfs1[,c("frequency","value")]
colnames(sfs1) <- c("frequency","count") # so it will work in the foldSFS function
sfs1_fold <- foldSFS(sfs1) # gets rid of the trailing zeroes (was folded already but wanted to elimiante extra zeroes)
sfs1_exclMono <- sfs1_fold[sfs1_fold$frequency!=0,]
sfs1_exclMono$label <- "1. 100% het filter (less stringent)"

sfs0.8$frequency <- as.numeric(sapply(strsplit(as.character(sfs0.8$variable),"_"),"[",2))
sfs0.8 <- sfs0.8[,c("frequency","value")]
colnames(sfs0.8) <- c("frequency","count") # so it will work in the foldSFS function
sfs0.8_fold <- foldSFS(sfs0.8) # gets rid of the trailing zeroes (was folded already but wanted to elimiante extra zeroes)
sfs0.8_exclMono <- sfs0.8_fold[sfs0.8_fold$frequency!=0,]
sfs0.8_exclMono$label <- "3. 80% het filter (more stringent)"

sfs0.7$frequency <- as.numeric(sapply(strsplit(as.character(sfs0.7$variable),"_"),"[",2))
sfs0.7 <- sfs0.7[,c("frequency","value")]
colnames(sfs0.7) <- c("frequency","count") # so it will work in the foldSFS function
sfs0.7_fold <- foldSFS(sfs0.7) # gets rid of the trailing zeroes (was folded already but wanted to elimiante extra zeroes)
sfs0.7_exclMono <- sfs0.7_fold[sfs0.7_fold$frequency!=0,]
sfs0.7_exclMono$label <- "4. 70% het filter (very stringent)"


sfs0.9$frequency <- as.numeric(sapply(strsplit(as.character(sfs0.9$variable),"_"),"[",2))
sfs0.9 <- sfs0.9[,c("frequency","value")]
colnames(sfs0.9) <- c("frequency","count") # so it will work in the foldSFS function
sfs0.9_fold <- foldSFS(sfs0.9) # gets rid of the trailing zeroes (was folded already but wanted to elimiante extra zeroes)
sfs0.9_exclMono <- sfs0.9_fold[sfs0.9_fold$frequency!=0,]
sfs0.9_exclMono$label <- "2. 90% het filter (less stringent)"


all <- rbind(sfs0.7_exclMono,sfs0.8_exclMono,sfs0.9_exclMono,sfs1_exclMono)
p <- ggplot(all,aes(x=frequency, y=count,fill=label))+
  geom_bar(stat="identity",position="dodge")+
  theme_bw()+
  scale_x_continuous(breaks=c(seq(1,length(sfs_exclMono$frequency))))+
  ggtitle("COM")
p
