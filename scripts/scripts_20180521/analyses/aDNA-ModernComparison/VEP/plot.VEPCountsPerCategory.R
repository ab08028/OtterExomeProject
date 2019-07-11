###### 
minInd=1
minDepth=1
################################ high cov ########################
HCdata.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/pseudoHaps/20190612-highcov-pseudoHaps/countsPerCategory/"
HCTotalCallable <- read.table(paste(HCdata.dir,"totalCallableCDSSitesPerInd.minDepth.",minDepth,".minInds.",minInd,".txt",sep=""),header=T)
HCSyn <- read.table(paste(HCdata.dir,"countsPerCategory.synonymous.minDepth.",minDepth,".minInds.",minInd,".txt",sep=""),header=T)
HCMis <- read.table(paste(HCdata.dir,"countsPerCategory.missense.minDepth.",minDepth,".minInds.",minInd,".txt",sep=""), header =T )
HCSG<- read.table(paste(HCdata.dir,"countsPerCategory.stopgained.minDepth.",minDepth,".minInds.",minInd,".txt",sep=""),header=T)
HCSG$state <- "highcov"
HCMis$state <- "highcov"
HCSyn$state <- "highcov"
HCTotalCallable$cat <- "allSites"
HCSG$cat <- "stopgained"
HCMis$cat <- "missense"
HCSyn$cat <- "synonymous"
HCSG_plusTotal = merge(HCTotalCallable[,c("sample","callableSites")],HCSG,by="sample",suffixes = c(".TotalSites",".Category"))
HCMis_plusTotal = merge(HCTotalCallable[,c("sample","callableSites")],HCMis,by="sample",suffixes = c(".TotalSites",".Category"))
HCSyn_plusTotal = merge(HCTotalCallable[,c("sample","callableSites")],HCSyn,by="sample",suffixes = c(".TotalSites",".Category"))

HCSG_plusTotal$FractionOfTotalCDS_TVOnly <- HCSG_plusTotal$nonRefCalls_TransversionsOnly / HCSG_plusTotal$callableSites.TotalSites

HCMis_plusTotal$FractionOfTotalCDS_TVOnly <- HCMis_plusTotal$nonRefCalls_TransversionsOnly / HCMis_plusTotal$callableSites.TotalSites

HCSyn_plusTotal$FractionOfTotalCDS_TVOnly <- HCSyn_plusTotal$nonRefCalls_TransversionsOnly / HCSyn_plusTotal$callableSites.TotalSites

allHC <- rbind(HCSG_plusTotal,HCMis_plusTotal,HCSyn_plusTotal)

allHC$group <- "NA"
allHC[grep("^A",allHC$sample),]$group <- "ancient"
allHC[grep("^A",allHC$sample,invert = T),]$group <- "modern"

################################ low cov ########################
LCdata.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/pseudoHaps/20190612-lowcov-pseudoHaps/countsPerCategory/"
LCTotalCallable <- read.table(paste(LCdata.dir,"totalCallableCDSSitesPerInd.minDepth.",minDepth,".minInds.",minInd,".txt",sep=""),header=T)
LCSyn <- read.table(paste(LCdata.dir,"countsPerCategory.synonymous.minDepth.",minDepth,".minInds.",minInd,".txt",sep=""),header=T)
LCMis <- read.table(paste(LCdata.dir,"countsPerCategory.missense.minDepth.",minDepth,".minInds.",minInd,".txt",sep=""), header =T )
LCSG<- read.table(paste(LCdata.dir,"countsPerCategory.stopgained.minDepth.",minDepth,".minInds.",minInd,".txt",sep=""),header=T)
LCSG$state <- "lowcov"
LCMis$state <- "lowcov"
LCSyn$state <- "lowcov"
LCTotalCallable$cat <- "allSites"
LCSG$cat <- "stopgained"
LCMis$cat <- "missense"
LCSyn$cat <- "synonymous"
LCSG_plusTotal = merge(LCTotalCallable[,c("sample","callableSites")],LCSG,by="sample",suffixes = c(".TotalSites",".Category"))
LCMis_plusTotal = merge(LCTotalCallable[,c("sample","callableSites")],LCMis,by="sample",suffixes = c(".TotalSites",".Category"))
LCSyn_plusTotal = merge(LCTotalCallable[,c("sample","callableSites")],LCSyn,by="sample",suffixes = c(".TotalSites",".Category"))

LCSG_plusTotal$FractionOfTotalCDS_TVOnly <- LCSG_plusTotal$nonRefCalls_TransversionsOnly / LCSG_plusTotal$callableSites.TotalSites

LCMis_plusTotal$FractionOfTotalCDS_TVOnly <- LCMis_plusTotal$nonRefCalls_TransversionsOnly / LCMis_plusTotal$callableSites.TotalSites

LCSyn_plusTotal$FractionOfTotalCDS_TVOnly <- LCSyn_plusTotal$nonRefCalls_TransversionsOnly / LCSyn_plusTotal$callableSites.TotalSites

allLC <- rbind(LCSG_plusTotal,LCMis_plusTotal,LCSyn_plusTotal)
allLC$group <- "NA"
allLC[grep("^A",allLC$sample),]$group <- "ancient"
allLC[grep("^A",allLC$sample,invert = T),]$group <- "modern"

############# plot ###########
p1 <- ggplot(allHC,aes(x=group,y=FractionOfTotalCDS_TVOnly,fill=cat))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~cat,scales="free")+
  theme_bw()+
  ggtitle(paste("Non-Downsampled Modern + Ancient\nFraction of Called Sites that are Derived Transversions (mapped to ferret)\nMinDepth=",minDepth," ; minInds=",minInd,sep=""))
p1
ggsave(paste(HCdata.dir,"highcov.fracOfCalledCDSSites.Derived.perCategory.pdf",sep=""),p1,height=5,width=7)
p2 <- ggplot(allLC,aes(x=group,y=FractionOfTotalCDS_TVOnly,fill=cat))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~cat,scales="free")+
  theme_bw()+
  ggtitle(paste("Downsampled Modern + Ancient\nFraction of Called Sites that are Derived Transversions (mapped to ferret)\nMinDepth=",minDepth," ; minInds=",minInd,sep=""))
p2
ggsave(paste(LCdata.dir,"lowcov.fracOfCalledCDSSites.Derived.perCategory.pdf",sep=""),p2,height=5,width=7)
