require(reshape2)
require(ggplot2)
require(RColorBrewer)
############# Plot **NEUTRAL **FSC format SFS in R ############
############### set up your colors -- keep this consistent across all plots ######
require(RColorBrewer)
colorPal=RColorBrewer::brewer.pal(n=6,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[7],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors

######### file specifics ######
genotypeDate="20181119"
dateOfProjection="20181221"
filter="hetFilter-0.75" # date you did the projection (so you can try different projections)
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/neutral/projection-",dateOfProjection,"-",filter,"/",sep="")
plot.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/SFS/",genotypeDate,"/easySFS_projection/projection-",dateOfProjection,"/",filter,"/neutral/",sep="")
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
  p <- ggplot(sfs_exclMono,aes(x=frequency, y=count,fill=pop))+
    geom_bar(stat="identity")+
    theme_bw()+
    scale_x_continuous(breaks=c(seq(1,length(sfs_exclMono$frequency))))+
    ggtitle(pop)+
    scale_fill_manual(values=colors[pop])
  ggsave(paste(plot.dir,pop,"folded.easySFS.proj.sfs.pdf",sep=""),device="pdf",height=5,width=7)
}

############# plot with transparent background ##########

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
  p <- ggplot(sfs_exclMono,aes(x=frequency, y=count,fill=pop))+
    geom_bar(stat="identity")+
    theme_classic()+
    scale_x_continuous(breaks=c(seq(1,length(sfs_exclMono$frequency))))+
    ggtitle(pop)+
    scale_fill_manual(values=colors[pop])+
    theme(panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA),text = element_text(size=20),axis.text =element_text(size=20))
  ggsave(paste(plot.dir,pop,"folded.easySFS.proj.sfs.transparent.pdf",sep=""),device="pdf",height=5,width=7)
}


############### 2D SFS ###################
############ need to update to automate ######## ## do this!
# so pop 0 is position 1 in pop Order (pythonic numbering). 

# test2D$pop1 <- as.character(row.names(test2D))
# sfs2D <- melt(test2D,id.vars = "pop1")
# head(sfs2D)
# colnames(sfs2D) <- c("pop1","pop0","count")
# sfs2D$frequency_pop0 <- as.numeric(sapply(strsplit(as.character(sfs2D$pop0),"_"),"[",2))
# sfs2D$frequency_pop1 <- as.numeric(sapply(strsplit(as.character(sfs2D$pop1),"_"),"[",2))
# head(sfs2D)
# ggplot(data = sfs2D, aes(x=frequency_pop0, y=frequency_pop1, fill=log10(count))) + 
#   geom_tile()+
#   scale_fill_gradientn(colours=rainbow(8))+
#   theme_bw()+
#   geom_text(aes(label = round(count,1))) +
#   scale_x_continuous(breaks=c(seq(0,length(sfs2D[sfs2D$frequency_pop1==0,]$frequency_pop0)-1)))+
#   scale_y_continuous(breaks=c(seq(0,length(sfs2D[sfs2D$frequency_pop0==0,]$frequency_pop1)-1)))

