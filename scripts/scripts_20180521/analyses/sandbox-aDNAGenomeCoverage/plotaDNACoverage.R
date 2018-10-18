# want to plot genome coverage
# but color regions red that are in capture, or something like that. 
# sum up coverage inside and outside capture regions? 
require(ggplot2)
require(reshape2)
require(dplyr)

date=20181015 # date of stats
todaysdate=format(Sys.Date(),format="%Y%m%d")
resultsDir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/"
plotoutDir=paste(resultsDir,"plots/mappingStats/aDNA/",date,"/",sep="") # where you want plots to do
tableoutDir=paste(resultsDir,"tables/mappingStats/aDNA/",date,"/",sep="")
dir.create(plotoutDir,showWarnings = F,recursive = T)
dir.create(tableoutDir,showWarnings = F,recursive = T)

################ PICARD HS METRICS ###################
picard <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/plmx_mappingStats/PicardHSMetrics/picard.HSMetrics.plmx.20181015.txt",header=T,sep="\t",stringsAsFactors = F)

View(picard)
# these are in reference to the sea otter genome

picard$category <- "single capture"
picard[grep("screen",picard$SAMPLE),]$category <- "uncaptured"
picard[grep("2CAP",picard$SAMPLE),]$category <- "double capture"

# definitions: https://broadinstitute.github.io/picard/picard-metric-definitions.html
variables <- c("TOTAL_READS","PF_BASES_ALIGNED",	"ON_BAIT_BASES",	"NEAR_BAIT_BASES",	"OFF_BAIT_BASES",	"ON_TARGET_BASES",	"PCT_SELECTED_BASES",	"PCT_OFF_BAIT",	"ON_BAIT_VS_SELECTED",	"MEAN_BAIT_COVERAGE",	"MEAN_TARGET_COVERAGE",	"PCT_USABLE_BASES_ON_TARGET",	"FOLD_ENRICHMENT",	"ZERO_CVG_TARGETS_PCT")
variables_sub <- c("TOTAL_READS","PCT_OFF_BAIT","MEAN_TARGET_COVERAGE","FOLD_ENRICHMENT",	"ZERO_CVG_TARGETS_PCT")
picard.melt <- melt(picard,id.vars = c("SAMPLE","category"),factorsAsStrings = T,measure.vars = variables)
picard.melt$reference <- "sea_otter_23May2016_bS9RH.deduped.99"


ggplot(picard.melt[picard.melt$variable %in% variables_sub,],aes(x=reorder(SAMPLE,-value),y=value,fill=category))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  coord_flip()+
  facet_wrap(~variable,scales = "free_x")+
  xlab("Library")+
  theme_bw()

########## Paleomix Mapping stats ############
plmx <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/plmx_mappingStats/20181015_allADNA/plmx.aDNA.summary.stats.allAncient.20181015.txt",header=T)

################### look at retained nts vs mapped nts #############
head(plmx)
colnames(plmx) <- c("SAMPLE","variable","reference","value")
plmx$category <- "single capture"
plmx[grep("screen",plmx$SAMPLE),]$category <- "uncaptured"
plmx[grep("2CAP",plmx$SAMPLE),]$category <- "double capture"
plmx$simpleRef <- NA
plmx[plmx$reference=="Mustela_putorius_furo.MusPutFur1.0.dna.toplevel",]$simpleRef <- "ferret"
plmx[plmx$reference=="sea_otter_23May2016_bS9RH.deduped.99",]$simpleRef <- "sea otter"
plmx[plmx$reference=="canFam3",]$simpleRef <- "dog"
plmx[plmx$reference=="ucsc.hg19",]$simpleRef <- "human"
# combine:
combo <- rbind.data.frame(plmx,picard.melt)

################ total fraction of sequenced bases that are aligned to sea otter ########
# this is the total sequenced bases passing filtering (unaligned)
nts <- combo[combo$variable=="seq_retained_nts",c("SAMPLE","value","category")]
# this is the total bases aligned to sea otter (note: bases not reads)
aligned <- combo[combo$variable=="PF_BASES_ALIGNED",c("SAMPLE","value","category")]
onBait <- combo[combo$variable=="ON_BAIT_BASES",c("SAMPLE","value","category")]

alnBase1 <- merge(nts, aligned,by=c("SAMPLE","category"),suffixes = c("",".aligned"))
alnBase2 <- merge(alnBase1,onBait,by=c("SAMPLE","category"),suffixes = c(".retained",".onBait"))
colnames(alnBase2) <- c("SAMPLE","category","retained_nts","aligned_nts","onBait_nts")
alnBase2$fracAlignedBases <- alnBase2$aligned_nts / alnBase2$retained_nts
alnBase2$fracOnBait <- alnBase2$onBait_nts / alnBase2$retained_nts

alnBase2.melt <- melt(alnBase2)
p0 <- ggplot(alnBase2.melt,aes(x=reorder(SAMPLE,-value),y=value,fill=category))+
  geom_bar(stat="identity")+
  coord_flip()+
  ggtitle("Comparison of capture performance when aligned to sea otter")+
  facet_wrap(~variable,scales="free_x")+
  ylab("")+
  xlab("")
p0
ggsave(paste(plotoutDir,"/ComparisonOfOnTargetRates.seaOtter.",todaysdate,".pdf",sep=""),p0,device="pdf",width = 11,height=8)
######## total fraction of sequenced bases that are on target #########
# this is the total sequenced bases passing filtering (unaligned)
nts <- combo[combo$variable=="seq_retained_nts",c("SAMPLE","value","category")]
# this is the total number of on-bait aligned bases 
onBait <- combo[combo$variable=="ON_BAIT_BASES",c("SAMPLE","value","category")]
# merge them: 
bases <- merge(nts,onBait,by=c("SAMPLE","category"),suffixes = c(".retained",".onBait"))
# get the on-bait fraction: 
bases$fracOnBait <- bases$value.onBait / bases$value.retained 
p1 <- ggplot(bases,aes(x=SAMPLE,y=fracOnBait,fill=category))+
  geom_bar(stat="identity")+
  coord_flip()+
  ylab("Fraction of total sequenced bases that are on-bait")
p1
ggsave(paste(plotoutDir,"/ComparisonOfOnTargetRates.seaOtter.",todaysdate,".pdf",sep=""),p0,device="pdf",width = 11,height=8)

