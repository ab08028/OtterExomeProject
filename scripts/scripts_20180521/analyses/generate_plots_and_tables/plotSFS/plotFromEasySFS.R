require(reshape2)
require(ggplot2)
require(RColorBrewer)
############# Plot FSC format SFS in R ############
test1D <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20180806/easySFS_projection/projection_20181121/fastsimcoal2/CA_MAFpop0.obs",skip = 1,header=T,stringsAsFactors = F) # this is a FOLDED sfs from EasySFS, but it has trailing zeroes
test2D <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20180806/easySFS_projection/projection_20181121/fastsimcoal2/neutral_jointMAFpop0_1.obs",skip = 1,header=T,stringsAsFactors = F)
  # how do i find out which are 0_1?

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
#################################### melt sfs ###################
sfs <- melt(test1D)
head(sfs)
# get frequency from the d0_0 
sfs$frequency <- as.numeric(sapply(strsplit(as.character(test1D_melt$variable),"_"),"[",2))
sfs <- sfs[,c("frequency","value")]
colnames(sfs) <- c("frequency","count") # so it will work in the foldSFS function
sfs_fold <- foldSFS(sfs) # gets rid of the trailing zeroes (was folded already but wanted to elimiante extra zeroes)
sfs_exclMono <- sfs_fold[sfs_fold$frequency!=0,]
ggplot(sfs_exclMono,aes(x=frequency, y=count))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_x_continuous(breaks=c(seq(1,length(sfs_exclMono$frequency))))

############### 2D SFS ###################
test2D$pop1 <- as.character(row.names(test2D))
sfs2D <- melt(test2D,id.vars = "pop1")
head(sfs2D)
colnames(sfs2D) <- c("pop1","pop0","count")
sfs2D$frequency_pop0 <- as.numeric(sapply(strsplit(as.character(sfs2D$pop0),"_"),"[",2))
sfs2D$frequency_pop1 <- as.numeric(sapply(strsplit(as.character(sfs2D$pop1),"_"),"[",2))
head(sfs2D)
ggplot(data = sfs2D, aes(x=frequency_pop0, y=frequency_pop1, fill=log10(count))) + 
  geom_tile()+
  scale_fill_gradientn(colours=rainbow(8))+
  theme_bw()+
  geom_text(aes(label = round(count,1))) +
  scale_x_continuous(breaks=c(seq(0,length(sfs2D[sfs2D$frequency_pop1==0,]$frequency_pop0)-1)))+
  scale_y_continuous(breaks=c(seq(0,length(sfs2D[sfs2D$frequency_pop0==0,]$frequency_pop1)-1)))
View(sfs2D)
