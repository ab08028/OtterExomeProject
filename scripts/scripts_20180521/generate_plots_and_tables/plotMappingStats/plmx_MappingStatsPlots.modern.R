require(ggplot2)
require(dplyr)
require(gridExtra)
require(scales)
############################## Modern DNA Mapping Assessment Plots ####################
#### This script explores the summary stats output by paleomix and preseq, etc.
#### where to save results:
resultsDir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/"
plotoutDir=paste(resultsDir,"plots/mappingStats/modernDNA/",sep="") # where you want plots to do
tableoutDir=paste(resultsDir,"tables/mappingStats/modernDNA/",sep="")
dir.create(plotoutDir,showWarnings = F,recursive = T)
dir.create(tableoutDir,showWarnings = F,recursive = T)

########################### 1. Paleomix Summary ###############################
# in bash, I selected the stats I wanted using sandbox--gatherstats.sh script
# These were then downloaded and can be input into R.

date=20180716 # date of stats
stats <- read.table(paste(resultsDir,"QC_Reports/plmx_mappingStats/",date,"/plmx.modernDNA.summary.stats.txt",sep=""),header=T)
dim(stats)
stats <- unique(stats)
dim(stats)
#head(stats)
# add short version of reference genome names: 
stats$simpleRef <- NA
stats[stats$reference=="Mustela_putorius_furo.MusPutFur1.0.dna.toplevel",]$simpleRef <- "ferret"
# for now, only mapping to mfur (but that could change)
#stats[stats$reference=="sea_otter_23May2016_bS9RH.deduped.99",]$simpleRef <- "sea otter"
#stats[stats$reference=="canFam3",]$simpleRef <- "dog"
#stats[stats$reference=="ucsc.hg19",]$simpleRef <- "human"
############## Add populations ##########
stats$population <- NA
stats[grep("_CA_",stats$sample),]$population <- "CA"
stats[grep("_AK_",stats$sample),]$population <- "AK"
stats[grep("_AL_",stats$sample),]$population <- "ALEUT"
stats[grep("_KUR_",stats$sample),]$population <- "KURIL"
stats[grep("_MED_",stats$sample),]$population <- "MEDNY"
stats[grep("_BER_",stats$sample),]$population <- "BERING"
######### Plot number of retained reads for each sample #######

# just plot seq retained reads:
p0 <- ggplot(data=stats[stats$statistic=="seq_retained_reads",],aes(x=reorder(sample,value),y=value,fill=population))+
  geom_bar(stat="identity")+
  coord_flip()+
  theme_bw()+
  ylab("Reads Passing Filters")+
  xlab("Sample")+
  facet_wrap(~statistic,scales="free_x")
p0
ggsave(paste(plotoutDir,"modernMappingStats.ReadsPerSample.",date,".pdf",sep=""),p0,device="pdf",width = 8,height=8)

# just plot seq retained reads:
p0.1 <- ggplot(data=stats[stats$statistic=="seq_retained_reads" & stats$value < 2e07,],aes(x=reorder(sample,value),y=value,fill=population))+
  geom_bar(stat="identity")+
  coord_flip()+
  theme_bw()+
  ylab("Reads Passing Filters")+
  xlab("Sample")+
  facet_wrap(~statistic,scales="free_x")+
  ggtitle("low coverage samples")
p0.1
ggsave(paste(plotoutDir,"modernMappingStats.ReadsPerSample.lowCoverageSamplesOnly.1e07.",date,".pdf",sep=""),p0.1,device="pdf",width = 8,height=8)

########### 2. Plot # of reads vs perc mapped #######
dfx <- stats[stats$statistic=="seq_retained_reads",c("sample","value","population")]
colnames(dfx) <- c("sample","Total_Reads","population")
dfy <- stats[stats$statistic=="hits_raw_frac",c("sample","value","population")]
colnames(dfy) <- c("sample","Frac_Mapped","population")
df <- merge(dfx,dfy,by=c("sample","population"))
dim(df)
p1 <- ggplot(df,aes(x=Total_Reads,y=Frac_Mapped,color=population))+
  geom_point()+
  theme_bw()
p1
ggsave(paste(plotoutDir,"modernMappingStats.TotalReads.vs.FracMapped.",date,".pdf",sep=""),p1,device="pdf",width = 8,height=8)
