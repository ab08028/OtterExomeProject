### 20200209: want to get total sites per individual with 1x coverage, 2x etc.
# then want to get total Tv SNPS per individual from later runs.

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

########## Get stats for SI:
# want to know: total reads covered by 1x in elut and mfur
# original only (not downsampled)
allDepths_melt[allDepths_melt$downsamp=="Original" & allDepths_melt$sampleID %in% c("A13_Elut_CA_AN_388_SN1_2CAP_screen","A29_Elut_CA_SM_30_SN2_CAP","A30_Elut_CA_SM_35_SN1_CAP"),]
# only aDNA samples:
# want 1-10 and then >10
# make your table:
tableForSI_1to10 <- allDepths_melt[allDepths_melt$downsamp=="Original" & allDepths_melt$sampleID %in% c("A13_Elut_CA_AN_388_SN1_2CAP_screen","A29_Elut_CA_SM_30_SN2_CAP","A30_Elut_CA_SM_35_SN1_CAP") & allDepths_melt$read_depth<=10 & allDepths_melt$read_depth>0,c("reference","sampleID","variable","value","read_depth")]


require(dplyr)    
tableForSI_greaterThan10 <- allDepths_melt[allDepths_melt$downsamp=="Original" & allDepths_melt$sampleID %in% c("A13_Elut_CA_AN_388_SN1_2CAP_screen","A29_Elut_CA_SM_30_SN2_CAP","A30_Elut_CA_SM_35_SN1_CAP") & 
allDepths_melt$read_depth>10,] %>%
  group_by(reference,sampleID)%>%
  summarise(value=sum(value),variable="SumGreaterThan10",read_depth=">10")

# combine:
tableForSICombo <- rbind(tableForSI_1to10,data.frame(tableForSI_greaterThan10))
write.table(tableForSICombo,paste(data.dir,"aDNA.ReadDepthSummaries.TotalAtEachReadDepth.NonInclusive.Upto10Thengt10.USEFORSI.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")

## get total reads with >1x
totalCoveredBases <- tableForSICombo %>%
  group_by(reference,sampleID)%>%
  summarise(totalCoveredBases = sum(value))
write.table(totalCoveredBases,paste(data.dir,"aDNA.TotalCoveredBases.forSI.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")

