require(ggplot2)
require(reshape2)
require(scales)
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/Depth/"
########### Plot depth ##################
refs=c("Mfur","Elut")
allDepths = data.frame()
for(ref in refs){
  bamList <- read.table(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/angsd.bamList.mappedto",ref,"fullpaths.txt",sep=""),stringsAsFactors = F)
  depthFilePerSample <- read.table(paste(data.dir,"angsd.Depths.mappedTo",ref,".depthSample",sep=""))
  colnames(depthFilePerSample) <- as.factor(paste("Depth_",seq(0,length(depthFilePerSample)-1),sep=""))
  depthFilePerSample$filename <- unlist(bamList)
  depthFilePerSample$baseFileName <- unlist(lapply(strsplit(depthFilePerSample$filename,"/"),tail,n=1))
  depthFilePerSample$sampleID <- unlist(lapply(strsplit(as.character(depthFilePerSample$baseFileName),"\\."),"[",1))
  depthFilePerSample$reference <- ref
  allDepths = rbind(allDepths,depthFilePerSample)
}
allDepths$downsamp <- "Original"
allDepths[grep("downsamp",allDepths$filename),]$downsamp <- "Downsampled"
allDepths_melt <- melt(allDepths,id.vars=c("reference","sampleID","downsamp"),measure.vars=paste("Depth_",seq(0,100),sep=""))
# want to exclude the 0 bin

allDepths_melt$read_depth <- as.numeric(unlist(lapply(strsplit(as.character(allDepths_melt$variable),"_"),"[",2)))
# separate by mfur / elut
# and per sample
for(ref in refs){
  for(sampleID in unique(allDepths_melt$sampleID)){
    p = ggplot(allDepths_melt[allDepths_melt$variable!="Depth_0" & allDepths_melt$sampleID==sampleID & allDepths_melt$reference==ref & allDepths_melt$value!=0,],aes(x=read_depth,y=value,fill=downsamp))+
      geom_bar(stat="identity",position="dodge")+
      ggtitle(paste(sampleID," mapped to ",ref,sep=""))+
      theme_bw()+
      theme(legend.title = element_blank(),legend.position = c(0.5,0.8),axis.title = element_text(size=14),legend.text = element_text(size=14))+
      scale_x_continuous(breaks = c(seq(0,100,by=5)))+
      scale_y_log10(labels=comma)+
      ylab("site count")+
      xlab("read depth")
      
    p
    ggsave(paste(data.dir,sampleID,".",ref,"depths.pdf",sep=""),p,height=5,width=5)
  }
}
