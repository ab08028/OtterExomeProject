require(ggplot2)
require(dplyr)
require(gridExtra)
require(scales)

##### Be careful -- these add stuff together if data points appear twice ###### 
############################## aDNA Assessment Plots ####################
#### This script explores the summary stats output by paleomix and preseq, etc.
#### where to save results:
date=20180730 # date of stats
todaysdate=format(Sys.Date(),format="%Y%m%d")
resultsDir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/"
plotoutDir=paste(resultsDir,"plots/mappingStats/aDNA/",date,"/",sep="") # where you want plots to do
tableoutDir=paste(resultsDir,"tables/mappingStats/aDNA/",date,"/",sep="")
dir.create(plotoutDir,showWarnings = F,recursive = T)
dir.create(tableoutDir,showWarnings = F,recursive = T)

############################# mt GENOME ##############
stats <- read.table(paste(resultsDir,"QC_Reports/plmx_mappingStats/",date,"/mtGenome/plmx.aDNA.summary.stats.mtGenome.incomplete.txt",sep=""),header=T,stringsAsFactors = F)
stats$simpleRef <- NA
stats[stats$reference=="mpt_ref_MusPutFur1.0_chrMT",]$simpleRef <- "ferret mtGenome"
stats[stats$reference=="elut_kenyoni_ref_ASM228890v2_chrMT",]$simpleRef <- "northern sea otter mtGenome"

############ set up df ready to plot ##############
statsToPlot1 <- stats[stats$statistic!="hits_raw" & stats$statistic!="hits_raw_frac",]
statsToPlot1$statistic <- factor(statsToPlot1$statistic)
levels(statsToPlot1$statistic) <- c("PCR Duplicate (frac)","Mapped Read Length\n(bp)", "Unique Mapped Reads\n(read count)","Unique Mapped Reads\n(frac)","Fraction of mtGenome\n>=1x","Fraction of mtGenome\n>=2x" ,"Fraction of mtGenome\n>=3x","Sequencing Reads")
statsToPlot1$plotOrder <- 1 # some lower number
statsToPlot1[statsToPlot1$sample %in% c("A9_Elut_CA_AN_388_SN1","A10_Elut_CA_NIC_1_SN1"),]$plotOrder <- 2
######### 0b Plot the Mb of genome coverage for samples that worked #########

samplesThatWorked <- statsToPlot1[statsToPlot1$statistic=="Fraction of mtGenome\n>=1x" & statsToPlot1$value>0,]$sample 
p0b <- ggplot(data=statsToPlot1[statsToPlot1$statistic %in% c("Fraction of mtGenome\n>=1x","Fraction of mtGenome\n>=2x","Fraction of mtGenome\n>=3x") & statsToPlot1$sample %in% samplesThatWorked,],aes(x=reorder(sample,plotOrder),y=value,fill=statistic))+
  geom_bar(stat="identity",position = "dodge")+
  ylab("Fraction of N. sea otter mtGenome Covered")+
  theme_bw()+
  coord_flip()+
  xlab("Sample")+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))
p0b
ggsave(paste(plotoutDir,"/Mb_NSeaOtter_MtGenome_MD_Coverage.onlySamplesThatWorked.",todaysdate,".pdf",sep=""),p0b,device="pdf",width = 8,height=5)


### Plot mapping stats  #####
p0c <- ggplot(data=statsToPlot1[statsToPlot1$statistic %in% c("Unique Mapped Reads\n(read count)","Unique Mapped Reads\n(frac)") & statsToPlot1$sample %in% samplesThatWorked,],aes(x=reorder(sample,plotOrder),y=value,fill=simpleRef))+
  geom_bar(stat="identity",position = "dodge")+
  ylab("")+
  theme_bw()+
  coord_flip()+
  xlab("Sample")+
  facet_wrap(~statistic,scales="free_x")+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14),legend.title = element_blank())
p0c
ggsave(paste(plotoutDir,"/Mb_NSeaOtter_Ferret_MtGenome_ReadsMapping.",todaysdate,".pdf",sep=""),p0c,device="pdf",width = 8,height=5)
