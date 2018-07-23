require(ggplot2)
require(dplyr)
require(gridExtra)
require(scales)

##### Be careful -- these add stuff together if data points appear twice ###### 
############################## aDNA Assessment Plots ####################
#### This script explores the summary stats output by paleomix and preseq, etc.
#### where to save results:
resultsDir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/"
plotoutDir=paste(resultsDir,"plots/mappingStats/aDNA/",sep="") # where you want plots to do
tableoutDir=paste(resultsDir,"tables/mappingStats/aDNA/",sep="")
dir.create(plotoutDir,showWarnings = F,recursive = T)
dir.create(tableoutDir,showWarnings = F,recursive = T)
########################### 1. Paleomix Summary ###############################
# in bash, I selected the stats I wanted using sandbox--gatherstats.sh script
# These were then downloaded and can be input into R.

date=20180722 # date of stats
stats <- read.table(paste(resultsDir,"QC_Reports/plmx_mappingStats/",date,"/plmx.aDNA.summary.stats.txt",sep=""),header=T,stringsAsFactors = F)
########## genome lengths for calculating 1x bases:
#mfur=2410758013
#cfam=2410976875
elut=2425597608 # only calculated the MD_001 for elut
#hg19=3137161264
########## CAUTION: if any duplicated entries, will stack up in bar plot -- so do unique always first #######
dim(stats)
stats <- unique(stats)
dim(stats)
# note how dim changes
# add short version of reference genome names: 
stats$simpleRef <- NA
stats[stats$reference=="Mustela_putorius_furo.MusPutFur1.0.dna.toplevel",]$simpleRef <- "ferret"
stats[stats$reference=="sea_otter_23May2016_bS9RH.deduped.99",]$simpleRef <- "sea otter"
stats[stats$reference=="canFam3",]$simpleRef <- "dog"
stats[stats$reference=="ucsc.hg19",]$simpleRef <- "human"

# make new rows with MB covered by 1,2,3x coverage statstics 
intMD <- stats[stats$statistic %in% c("MD_001","MD_002","MD_003"),]
intMD$value <- (intMD$value * elut) / 1e6
intMD[intMD$statistic=="MD_001",]$statistic <- "Mb >= 1x" # note this is mb of elut genome
intMD[intMD$statistic=="MD_002",]$statistic <- "Mb >= 2x"
intMD[intMD$statistic=="MD_003",]$statistic <- "Mb >= 3x"

# and add it in:
stats <- rbind(stats,intMD)
tail(stats)
######### 0. Plot number of retained reads for each sample #######

# just plot seq retained reads:
p0 <- ggplot(data=stats[stats$statistic=="seq_retained_reads",],aes(x=reorder(sample,value),y=value))+
  geom_bar(stat="identity")+
  coord_flip()+
  theme_bw()+
  ylab("Reads Passing Filters")+
  xlab("Sample")+
  facet_wrap(~statistic,scales="free_x")
p0

######### 0b Plot the Mb of genome coverage at 1/2x per sample #########
p0b <- ggplot(data=stats[stats$statistic %in% c("Mb >= 1x"),],aes(x=sample,y=value,fill=statistic))+
  geom_bar(stat="identity",position = "dodge")+
  ylab("Mb of sea otter genome with >= 1x coverage")+
  theme_bw()
p0b
ggsave(paste(plotoutDir,"MbSeaOtterGenome_1xCoverage.",date,".pdf",sep=""),p0b,device="pdf",width = 8,height=8)



######### Function: Plot stats by sample/reference ############
############ *** CAUTION **** If there are any duplicated entries in the dataframe, they will stack on top of each other ; added border (color=i) to check for this *be careful in all plots* #####

plotStatsFunc <- function(x, na.rm = TRUE, ...) {
  p <- list()
  nm <- unique(x$sample)
  for (i in seq_along(nm)) {
    p[[i]] <- ggplot(data=unique(filter(x,sample==nm[i])),aes(x=simpleRef,y=value,fill=simpleRef,color=i))+
          geom_bar(stat="identity") +
          theme_bw()+
          xlab("")+
          ylab("")+
          theme(legend.position = "none")+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          facet_wrap(~statistic,scales="free",nrow=1) +
          ggtitle(nm[i]) }
  return(p)
}

## Pick stats I want to plot with this
statsToPlot <- stats[stats$statistic!="seq_retained_reads" & stats$statistic!="hits_raw" & stats$statistic!="hits_raw_frac"  & !(stats$statistic %in% c("MD_001","MD_002","MD_003")) & !(stats$statistic %in% c("Mb >= 1x","Mb >= 2x","Mb >= 3x")),]

## drop old levels and rename new levels to be more interpretable:
statsToPlot$statistic <- factor(statsToPlot$statistic)
levels(statsToPlot$statistic) <- c("PCR Duplicates","Mapped Read Length\n(bp)", "Unique Mapped Reads\n(read count)","Unique Mapped Reads\n(frac)")

## Make plots for each sample
p1 <- plotStatsFunc(statsToPlot)
# put in a grid (one sample per row)
p1_grid <- do.call(grid.arrange,p1) # # do call calls grid.arrange for everything in p https://stackoverflow.com/questions/9315611/grid-of-multiple-ggplot2-plots-which-have-been-made-in-a-for-loop
# how do I save these?
ggsave(paste(plotoutDir,"aDNAMappingStats",date,".pdf",sep=""),p1_grid,device="pdf",width = 8,height=8)


## make plots for each statistic mapped to ferret

p1b <- ggplot(data=statsToPlot[statsToPlot$simpleRef=="sea otter",],aes(x=reorder(sample,value),y=value))+
  geom_bar(stat="identity")+
  facet_wrap(~statistic,scales = "free_x")+
  coord_flip()+
  xlab("")+
  ylab("")+
  ggtitle("aDNA Mapping Stats")+
  theme_bw()
p1b
ggsave(paste(plotoutDir,"aDNA.AllSamples.Stats.MfurRef.",date,".pdf",sep=""),p1b,device="pdf",width = 8,height=8)

########### 2. Plot # of reads vs perc mapped #######
dfx <- stats[stats$statistic=="seq_retained_reads",c("sample","value")]
colnames(dfx) <- c("sample","Total_Reads")
dfy <- stats[stats$statistic=="hits_raw_frac",c("sample","value","simpleRef")]
colnames(dfy) <- c("sample","Frac_Mapped","simpleRef")
df <- merge(dfx,dfy,by="sample")
dim(df)
p2 <- ggplot(df,aes(x=Total_Reads,y=Frac_Mapped,group=interaction(simpleRef,sample),color=sample))+
  geom_point()+
  facet_wrap(~simpleRef)
p2
ggsave(paste(plotoutDir,"aDNA_ReadCount-vs-PercMapped.",date,".pdf",sep=""),p2,device="pdf",width = 8,height=8)

########### 3. Compare my results to Avila-Arcos study #######
arcos <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/ComparisonToArcos_aDNA/Arcos_Table1_Table2_aDNACaptureStats.txt",header=T)
colnames(arcos) <- c("sample","age","seq_retained_reads","hits_unique","hits_unique_frac","technique","pcr_duplicates_perc","Mb >= 1x")
# pcr dups is not the same as clonailty of hits; don't compare 
# convert sequencing reads to millions
arcos$seq_retained_reads <- arcos$seq_retained_reads * 1e6
arcos$study <- "Avila-Arcos"
arcos.melt <- melt(arcos)
arcos.melt.comp <- arcos.melt[,c("sample","variable","value","study","technique")]
colnames(arcos.melt.comp) <- c("sample","statistic","value","study","technique")
# exclude age:
arcos.melt.comp <- arcos.melt.comp[arcos.melt.comp$statistic!="age" & arcos.melt.comp$statistic!="pcr_duplicates_perc",]

# compare to sea otter mapped seqs and sequenced reads
stat.comp <- stats[stats$reference=="sea_otter_23May2016_bS9RH.deduped.99" | stats$reference=="noREF",c("sample","statistic","value")]

# exclude some stats:
stat.comp <- stat.comp[stat.comp$statistic %in% c("seq_retained_reads","hits_unique", "hits_unique_frac","Mb >= 1x"),]
stat.comp$study <- "this study"
stat.comp$technique <- "Sea Otter Capture"

compElutArcos <- rbind(arcos.melt.comp,stat.comp)
# for now just compare mybaits to sea otter capture 
p3 <- ggplot(compElutArcos[compElutArcos$technique!="1_WISC",],aes(x=sample,y=value,fill=technique))+
  geom_bar(stat="identity",position = "dodge")+
  facet_wrap(~statistic,scales="free_x")+
  coord_flip()+
  ggtitle("Comparing to Avila-Arcos et al.")
p3
ggsave(paste(plotoutDir,"aDNA_compareToAvilaArcosStudy.",date,".pdf",sep=""),p3,device="pdf",width = 8,height=8)

########################### X. Preseq Library Complexity ###############################
# run preseq on the histograms output by paleomix (have to specify that you want them)
# put samples together? differs by reference