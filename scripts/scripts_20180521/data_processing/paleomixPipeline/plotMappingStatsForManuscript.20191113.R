require(RColorBrewer)
colorPal=RColorBrewer::brewer.pal(n=8,name = "Dark2")
# make Baja brown
colors=list(CA=colorPal[1],BAJ=colorPal[7],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors

##### want to make some SI sequencing plots #####
plmx <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/plmx_mappingStats/20180730/plmx.modernDNA.summary.stats.txt",header=T)
plmxBaja <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/plmx_mappingStats/20181103_BajaPlusDogs/plmx.modernDNA.summary.stats.BajaPlusDogs.txt",header=T)
plmxCombo <- rbind(plmx,plmxBaja)

# [1] "hits_clonality"     "hits_length"        "hits_raw"           "hits_raw_frac"     
# [5] "hits_unique"        "hits_unique_frac"   "seq_retained_reads"
# want to make some plots
stats=c("hits_unique","hits_clonality","seq_retained_reads")

# don't plot the RWAB pilot samples since excluded earlier on
plmxCombo_noRWAB <- plmxCombo[grep("RWAB",plmxCombo$sample,invert = T),]
plmxCombo_noRWAB_noCfam <- plmxCombo_noRWAB[grep("Cfam",plmxCombo_noRWAB$sample,invert=T),]
require(ggplot2)
plmxCombo_noRWAB_noCfam_selectStats <- plmxCombo_noRWAB_noCfam[plmxCombo_noRWAB_noCfam$statistic %in% stats,]

plmxCombo_noRWAB_noCfam_selectStats$statistic <- factor(plmxCombo_noRWAB_noCfam_selectStats$statistic,levels=c("seq_retained_reads","hits_clonality","hits_unique"))
###### add population info:
plmxCombo_noRWAB_noCfam_selectStats$population <- unlist(lapply(strsplit(as.character(plmxCombo_noRWAB_noCfam_selectStats$sample),"_"),"[",3))
### combine BER and MED into one pop
plmxCombo_noRWAB_noCfam_selectStats$population2 <- plmxCombo_noRWAB_noCfam_selectStats$population
plmxCombo_noRWAB_noCfam_selectStats[plmxCombo_noRWAB_noCfam_selectStats$population %in% c("BER","MED"),]$population2 <- "COM"

### want to rename the statistics
plmxCombo_noRWAB_noCfam_selectStats$statLabel <- NA
plmxCombo_noRWAB_noCfam_selectStats[plmxCombo_noRWAB_noCfam_selectStats$statistic=="seq_retained_reads",]$statLabel <- "total seq. reads"
plmxCombo_noRWAB_noCfam_selectStats[plmxCombo_noRWAB_noCfam_selectStats$statistic=="hits_clonality",]$statLabel <- "clonality"
plmxCombo_noRWAB_noCfam_selectStats[plmxCombo_noRWAB_noCfam_selectStats$statistic=="hits_unique",]$statLabel <- "uniq. mapped reads"
plmxCombo_noRWAB_noCfam_selectStats$statLabel <- factor(plmxCombo_noRWAB_noCfam_selectStats$statLabel,levels=c("total seq. reads","clonality","uniq. mapped reads"))
p1 <- ggplot(plmxCombo_noRWAB_noCfam_selectStats,aes(x=reorder(sample,value),y=value,fill=population2))+
  geom_col()+
  coord_flip()+
  theme_bw()+
  facet_wrap(~statLabel,ncol = 3,scales="free_x")+
  xlab("")+
  ylab("")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c(colors['CA'],colors['BAJ'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))
p1

todaysdate=format(Sys.Date(),format="%Y%m%d")
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/mappingStats/modernDNA/sequencingReads.clonality.hitstoMfur.AllModern.IncludesBaja.",todaysdate,".pdf",sep=""),p1,height=15,width=10)
####### color by whether sample was filtered due to missing GTs down stream ######
filtered <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/plmx_mappingStats/20180730/samplesRemovedDueToHighMissingGTslater.txt",header=F)

plmxCombo_noRWAB_noCfam_selectStats$label <- "not filtered due to missing genotypes"
plmxCombo_noRWAB_noCfam_selectStats[plmxCombo_noRWAB_noCfam_selectStats$sample %in% filtered$V1,]$label <- "filtered due to missing genotypes"



p2 <- ggplot(plmxCombo_noRWAB_noCfam_selectStats,aes(x=reorder(sample,value),y=value,fill=label))+
  geom_col()+
  coord_flip()+
  theme_bw()+
  facet_wrap(~statistic,ncol = 3,scales="free_x")+
  xlab("")+
  ylab("")
p2
