require(ggplot2)
require(scales)
######### aDNA plotting stats for supplement: A13, A29, A30 only:
# do we need to plot all the failed ones as well? Not sure how that all goes. # can add if asked?

plmxStats <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/plmx_mappingStats/20190926_fornewplots_aDNA.plusA29.A30/allaDNAPlmxStatsTogether.plusA29.A30.txt",header=T,stringsAsFactors = F)
plot.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/A13_A29_A30_StatReportForManuscript/"
head(plmxStats)
plmxStats$filter <- "not-used"
# the ones chosen to go forward are A13,A29 and A30
plmxStats[plmxStats$sampleShortName %in% c("A13","A29","A30"),]$filter <- "used"
plmxStats$sampleShortName <- unlist(lapply(strsplit(plmxStats$sample,"_"),"[",1))
plmxStats$descriptive <- unlist(lapply(strsplit(plmxStats$sample,"_Elut_"),"[",2))
# got to add lib blank:
plmxStats[grepl("Blank",plmxStats$sample),]$descriptive <- "libraryBlank"
######### NOTE there are multiple measures PER sample because of different references 
# just use ferret reference? or sea otter?

####### first just want to show why I chose A13,A29 and A30: #######
p1 <- ggplot(plmxStats[plmxStats$statistic=="hits_unique" & plmxStats$reference=="sea_otter_23May2016_bS9RH.deduped.99",],aes(x=sample,y=value,fill=filter))+
  geom_bar(stat="identity",position="dodge")+
  theme_bw()+
  xlab("")+
  ylab("Number of reads that map uniquely to sea otter genome")+
  coord_flip()+
  ggtitle("The three ancient samples chosen for futher analysis (blue)\nhave the most reads mapping to the sea otter genome")+
  theme(legend.position = "none")+
  scale_y_continuous(labels=comma)
p1
ggsave(paste(plot.dir,"uniqueReadCount.allaDNASamples.pdf",sep=""),height=5,width=9)
######
######## show more properties of those good 3: ########
# note A29 was first screened as A19 
# but I don't have canfam or human info for A29 (was too slow), but showed initial low levels of contam. But maybe I don't need to do show this. (just have a backup)
statsToPlot1 <- c("hits_unique_frac","hits_unique")
p2Temp <- ggplot(plmxStats[plmxStats$statistic %in% statsToPlot1 & plmxStats$filter=="used",],aes(x=sampleShortName,y=value,fill=reference))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~statistic,scales="free")+
  coord_flip()+
  theme_bw()
p2Temp
ggsave(paste(plot.dir,"notNeeded.MappingToDogHuman.A30.A13.MissingA29ButCouldUseA19instead.pdf",sep=""),p2Temp,height=5,width=9)

########## Show read distribution for A13,A29,A30 #########
statsToPlot2 <- c("hits_unique","hits_clonality","hits_length")
refs=c("sea_otter_23May2016_bS9RH.deduped.99","Mustela_putorius_furo.MusPutFur1.0.dna.toplevel")
### only plot top 3:
p3a <- ggplot(plmxStats[plmxStats$reference %in% refs & plmxStats$statistic %in% statsToPlot2 & plmxStats$filter=="used",],aes(x=sampleShortName,y=value,fill=reference))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~statistic,scales="free")+
  theme_bw()+
  xlab("")+
  ylab("")
p3a
ggsave(paste(plot.dir,"A13.A29.A30.clonality.hitsLengths.TotalLength.ElutVsMfur.pdf",sep=""),p3a,height=5,width=9)

# just showing diff in mapping between ancient/modern: (but note in a plot made elsewhere when you add in modern samples they DON'T show the differnce-- PHEW PHEW PHEW.)
statsToPlot3 <- c("hits_unique")
refs=c("sea_otter_23May2016_bS9RH.deduped.99","Mustela_putorius_furo.MusPutFur1.0.dna.toplevel") 
p3b <- ggplot(plmxStats[plmxStats$reference %in% refs & plmxStats$statistic %in% statsToPlot3 & plmxStats$filter=="used",],aes(x=sampleShortName,y=value,fill=reference))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~statistic,scales="free")+
  theme_bw()+
  xlab("")+
  ylab("reads")+
  scale_y_continuous(labels=comma)
p3b
ggsave(paste(plot.dir,"A13.A29.A30.ElutVsMfurHits.pdf",sep=""),p3b,height=4,width=7)
