require(ggplot2)
genotypeDate="20181119"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/compareMissense_Syn/",genotypeDate,"/countsOfGenotypesPerIndividual/",sep="")
colorPal=RColorBrewer::brewer.pal(n=8,name = "Dark2")
# make Baja brown
colors=list(CA=colorPal[1],BAJ=colorPal[7],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors
pops=c("CA","AK","AL","COM","KUR") # excluding BAJA because too low sample size
date="20191014" # annoying -- maybe don't output this in future
categories=c("syn","missense")
missingness=c("0.8","0.90","0.95","1") # this is max missingness threshold from vcftools which is unhelpfully named. it's really the min call rate, so 0 refers to any amount of missingness being allowed  (min call rate is 0); 0.8 allows 20% missingness; 1 allows no missingness. very confusing nomenclature, I know. Leaving 0 out because already plotted using old script and it has a slightly different filenames structure -- look in "olderScripts" to see plottin g of missingness  = 0 (all missing allowed)
allCountsSynMis <- data.frame()
for(missing in missingness){
  for(cat in categories){
    counts <- read.table(paste(data.dir,"/minCallRate_",missing,"/",cat,".countsPerIndividual.minCall.",missing,".countOfHomAltRefHet.txt",sep=""),header=T) ### DONT USE! this has no missingness! use missingness filters (see below for example)
    counts$derivedAlleles <- (2*counts$HomAltCount) + counts$HetCount # note -- Hom ref counts are very low because these are all snps, many of which are fixed 1/1 against ferret (elevates count of hom alt a lot); monomorphic ref sites will contain all the 0/0s.
    counts$category <- cat
    counts$missingLevel <- missing
    # normalize by missingness per individual
    # doing the list.files thing to avoid needin the specific date on the file
    totalCallableCDSPerInd <- read.table(list.files(path=paste(data.dir,"/callableCDSPerIndividual/","minCallRate_",missing,sep=""),pattern=paste("callableCDSSitesPerIndividual.*.txt",sep=""),full.names = T),header=T)
    counts_PlusCallable <- merge(counts,totalCallableCDSPerInd,by.x = "individual",by.y="id")
    avgCallableCDS = mean(totalCallableCDSPerInd$CalledCount)
    counts_PlusCallable$derivedAlleles = (2*counts$HomAltCount) + counts$HetCount
    counts_PlusCallable$derivedAlleles_rescaled <- (counts_PlusCallable$derivedAlleles/counts_PlusCallable$CalledCount)*avgCallableCDS
    counts_PlusCallable$homAlt_rescaled <- (counts_PlusCallable$HomAltCount / counts_PlusCallable$CalledCount)*avgCallableCDS
    counts_PlusCallable$het_rescaled <- (counts_PlusCallable$HetCount / counts_PlusCallable$CalledCount)*avgCallableCDS
    allCountsSynMis = rbind(allCountsSynMis,counts_PlusCallable)
    
  }
}


##### get populations: ####
allCountsSynMis$pop1 <- unlist(lapply(strsplit(as.character(allCountsSynMis$individual),"_"),"[",3))
# make Commanders one population :
allCountsSynMis$pop2 <- allCountsSynMis$pop1
allCountsSynMis[allCountsSynMis$pop1 %in% c("BER","MED"),]$pop2 <- "COM"

######## order factors ######

allCountsSynMis$category <- factor(allCountsSynMis$category,levels=c("syn","missense"))
allCountsSynMis$pop2 <- factor(allCountsSynMis$pop2,levels=c("CA","BAJ","AK","AL","COM","KUR"))

######### plot with facets now ########
p1 <- ggplot(allCountsSynMis,aes(x=pop2,y=derivedAlleles_rescaled,fill=pop2))+
  geom_boxplot()+
  geom_point(alpha=0.4)+
  facet_wrap(~missingLevel~category,scales="free",ncol=2)+
  theme_bw()+
  scale_fill_manual(values=c(colors['CA'],colors['BAJ'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme(legend.position = "None")+
  ggtitle("derived alleles (normalized)")+
  xlab("")
p1
ggsave(paste(data.dir,"/derivedAlleles.Perpop.MultMissingFilts.withBaja.pdf",sep=""),p1,height=15,width=8)
################# drop baja ######################
p2 <- ggplot(allCountsSynMis[allCountsSynMis$pop2!="BAJ",],aes(x=pop2,y=derivedAlleles_rescaled,fill=pop2))+
  geom_boxplot()+
  geom_point(alpha=0.4)+
  facet_wrap(~missingLevel~category,scales="free",ncol=2)+
  theme_bw()+
  scale_fill_manual(values=c(colors['CA'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme(legend.position = "None")+
  ggtitle("derived alleles (normalized)")+
  theme(text=element_text(size=14))+
  xlab("")
p2
ggsave(paste(data.dir,"/derivedAlleles.Perpop.MultMissingFilts.noBaja.pdf",sep=""),p2,height=15,width=8)

############### Plot hom - alt and hets ##############
p3 <- ggplot(allCountsSynMis[allCountsSynMis$pop2!="BAJ",],aes(x=pop2,y=homAlt_rescaled,fill=pop2))+
  geom_boxplot()+
  geom_point(alpha=0.4)+
  facet_wrap(~missingLevel~category,scales="free",ncol=2)+
  theme_bw()+
  scale_fill_manual(values=c(colors['CA'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme(legend.position = "None")+
  ggtitle("hom-alt genotypes (normalized)")+
  theme(text=element_text(size=14))+
  xlab("")
p3
ggsave(paste(data.dir,"/homAlt.Perpop.MultMissingFilts.noBaja.pdf",sep=""),p3,height=15,width=8)
###### with baja
p4 <- ggplot(allCountsSynMis,aes(x=pop2,y=homAlt_rescaled,fill=pop2))+
  geom_boxplot()+
  geom_point(alpha=0.4)+
  facet_wrap(~missingLevel~category,scales="free",ncol=2)+
  theme_bw()+
  scale_fill_manual(values=c(colors['CA'],colors['BAJ'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme(legend.position = "None")+
  ggtitle("hom-alt genotypes (normalized)")+
  theme(text=element_text(size=14))+
  xlab("")
p4
ggsave(paste(data.dir,"/homAlt.Perpop.MultMissingFilts.pdf",sep=""),p4,height=15,width=8)

#### plot hets, no baja:
p4b <- ggplot(allCountsSynMis[allCountsSynMis$pop2!="BAJ",],aes(x=pop2,y=het_rescaled,fill=pop2))+
  geom_boxplot()+
  geom_point(alpha=0.4)+
  facet_wrap(~missingLevel~category,scales="free",ncol=2)+
  theme_bw()+
  scale_fill_manual(values=c(colors['CA'],colors['BAJ'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme(legend.position = "None")+
  ggtitle("het genotypes (normalized)")+
  theme(text=element_text(size=14))+
  xlab("")
p4b
ggsave(paste(data.dir,"/hets.Perpop.MultMissingFilts.noBaja.pdf",sep=""),p4b,height=15,width=8)

################### try getting rid of extreme outliers #############
head(allCountsSynMis)
# look at dist of called:
# just total called CDS:
calledCDSOnly <- unique(allCountsSynMis[,c("individual","CalledCount","missingLevel","pop2")])
ggplot(calledCDSOnly,aes(x=individual,y=CalledCount))+
  geom_point()+
  facet_wrap(~missingLevel,ncol=1,scales="free")
require(dplyr)
calledCDSOnly <- calledCDSOnly %>% 
  group_by(missingLevel) %>%
  mutate(meanCallable=mean(CalledCount),sdCallable=sd(CalledCount)) # adding in mean and sd callable sites for each group of missingness 
# label those that are < 1sd from mean callable:
calledCDSOnly$label <- "PASS"
calledCDSOnly[calledCDSOnly$CalledCount < (calledCDSOnly$meanCallable - calledCDSOnly$sdCallable), ]$label <- "FAIL"
pCallable <- ggplot(calledCDSOnly,aes(x=individual,y=CalledCount,color=label))+
  geom_point()+
  facet_wrap(~missingLevel,ncol=1,scales="free")+
  ggtitle("Individuals with < 1 sd of mean callable sites (for each missingness thresholdhold) in red")+
  theme_bw()
ggsave(paste(data.dir,"/CallableCDSSitesPerIndividual.Below1SDInRed.pdf",sep=""),pCallable,height=5,width=8)

# get the names of those individuals and then exclude them from next plot:
failingInds <- calledCDSOnly[calledCDSOnly$label=="FAIL",]$individual


p2_alt <- ggplot(allCountsSynMis[allCountsSynMis$pop2!="BAJ" & !(allCountsSynMis$individual %in% failingInds),],aes(x=pop2,y=derivedAlleles_rescaled,fill=pop2))+
  geom_boxplot()+
  geom_point(alpha=0.4)+
  facet_wrap(~missingLevel~category,scales="free",ncol=2)+
  theme_bw()+
  scale_fill_manual(values=c(colors['CA'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme(legend.position = "None")+
  ggtitle("derived alleles (normalized)")+
  theme(text=element_text(size=14))+
  xlab("")
p2_alt
ggsave(paste(data.dir,"/derivedAlleles.Perpop.MultMissingFilts.noBaja.noOutliers.pdf",sep=""),p2_alt,height=15,width=8)

p3_alt <- ggplot(allCountsSynMis[allCountsSynMis$pop2!="BAJ" & !(allCountsSynMis$individual %in% failingInds),],aes(x=pop2,y=homAlt_rescaled,fill=pop2))+
  geom_boxplot()+
  geom_point(alpha=0.4)+
  facet_wrap(~missingLevel~category,scales="free",ncol=2)+
  theme_bw()+
  scale_fill_manual(values=c(colors['CA'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme(legend.position = "None")+
  ggtitle("hom-alt genotypes (normalized)")+
  theme(text=element_text(size=14))+
  xlab("")
p3_alt
ggsave(paste(data.dir,"/homAlt.Perpop.MultMissingFilts.noBaja.NoOutliers.pdf",sep=""),p3_alt,height=15,width=8)
############ Plot called vs rescaled #########
p5a <- ggplot(allCountsSynMis[allCountsSynMis$pop2!="BAJ" & allCountsSynMis$missingLevel!=1,],aes(x=as.numeric(CalledCount),y=derivedAlleles_rescaled,color=pop2))+
  geom_point()+
  facet_wrap(~missingLevel~category,scales="free")+
  theme_bw()+
  scale_color_manual(values=c(colors['CA'],colors['BAJ'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme(legend.position = "None")+
  ggtitle("derived alleles vs called cds sites ")+
  theme(text=element_text(size=14))+
  xlab("")
p5a
ggsave(paste(data.dir,"/derivedAlleles.vs.CalledSites.multfilt.pdf",sep=""),p5a,height=5,width=7)


p5b <- ggplot(allCountsSynMis[allCountsSynMis$pop2!="BAJ" & allCountsSynMis$missingLevel==0.8,],aes(x=as.numeric(CalledCount),y=homAlt_rescaled,color=pop2))+
  geom_point()+
  facet_wrap(~missingLevel~category,scales="free")+
  theme_bw()+
  scale_color_manual(values=c(colors['CA'],colors['BAJ'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme(legend.position = "None")+
  ggtitle("hom-alt GTs vs called cds sites ")+
  theme(text=element_text(size=14))+
  xlab("")
p5b
ggsave(paste(data.dir,"/homAltRescaled.vs.CalledSites.pdf",sep=""),p5b,height=5,width=7)


p5c <- ggplot(allCountsSynMis[allCountsSynMis$pop2!="BAJ" & allCountsSynMis$missingLevel==0.8,],aes(x=as.numeric(CalledCount),y=het_rescaled,color=pop2))+
  geom_point()+
  facet_wrap(~missingLevel~category,scales="free")+
  theme_bw()+
  scale_color_manual(values=c(colors['CA'],colors['BAJ'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme(legend.position = "None")+
  ggtitle("het GTs vs called cds sites ")+
  theme(text=element_text(size=14))+
  xlab("")
p5c
ggsave(paste(data.dir,"/hetRescaled.vs.CalledSites.pdf",sep=""),p5c,height=5,width=7)

############### Playing with beta #########################

beta = 
############## SANDBOX: playing with linear models ############
# http://www.sthda.com/english/articles/40-regression-analysis/168-multiple-linear-regression-in-r/
################## model1b : derived alleles frac vs called count ##########
testDF = allCountsSynMis[allCountsSynMis$category=="syn" & allCountsSynMis$missingLevel==0.8 & allCountsSynMis$pop2!="BAJ",]
testDF$derivedAlleles_frac <- testDF$derivedAlleles/testDF$CalledCount
# want to get lm of that
model1b <- lm(derivedAlleles_frac~CalledCount,data=testDF)
summary(model1b)
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 2.137e-02  3.389e-04   63.05   <2e-16 ***
#CalledCount 2.113e-10  1.746e-11   12.10   <2e-16 ***
ggplot(testDF,aes(x=CalledCount,y=derivedAlleles_frac))+
  geom_point()+
  geom_smooth(method="lm")
summary(model1b)
testDF$test1 <- testDF$derivedAlleles_frac - (2.113e-10*testDF$CalledCount) # aha so subtract off the influence of added called sites! 
ggplot(testDF,aes(x=CalledCount,y=test1))+
  geom_point()+
  geom_smooth(method="lm")+
  ggtitle("subtract off mx term")
meanCalled=mean(testDF$CalledCount)

testDF$test1_scaled <- testDF$test1 * meanCalled
ggplot(testDF,aes(x=CalledCount,y=test1_scaled))+
  geom_point()+
  geom_smooth(method="lm")+
  ggtitle("subtract off mx term")

ggplot(testDF,aes(x=pop2,y=test1_scaled))+
  geom_boxplot()
############## how does this affect hom alt? ##########


lm(homAlt_rescaled~CalledCount,data=allCountsSynMis[allCountsSynMis$category=="missense" & allCountsSynMis$missingLevel==0.8,])
#(Intercept)  CalledCount  
#9.670e+04    4.006e-04
lm(homAlt_rescaled~CalledCount,data=allCountsSynMis[allCountsSynMis$category=="syn" & allCountsSynMis$missingLevel==0.8,])
# (Intercept)  CalledCount  
# 2.077e+05    1.995e-03  

##### try excluding ####
p5 <- ggplot(allCountsSynMis[allCountsSynMis$pop2!="BAJ" & allCountsSynMis$missingLevel==0.8,],aes(x=as.numeric(CalledCount),y=derivedAlleles_rescaled,color=pop2))+
  geom_point()+
  facet_wrap(~missingLevel~category,scales="free")+
  theme_bw()+
  scale_color_manual(values=c(colors['CA'],colors['BAJ'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme(legend.position = "None")+
  ggtitle("derived alleles vs called cds sites ")+
  theme(text=element_text(size=14))+
  xlab("")
p5

### plot called count:
ggplot(totalCallableCDSPerInd,aes(x=id,y=CalledCount))+
  geom_point()
sd(totalCallableCDSPerInd$CalledCount)
mean(totalCallableCDSPerInd$CalledCount)
# get inds that are too low: 
goodInds <- totalCallableCDSPerInd[totalCallableCDSPerInd$CalledCount>1.92e7,]
badInds <- totalCallableCDSPerInd[totalCallableCDSPerInd$CalledCount<=1.92e7,]
ggplot(allCountsSynMis[allCountsSynMis$individual %in% goodInds$id & allCountsSynMis$missingLevel==0.8 & allCountsSynMis$pop2!="BAJ",],aes(x=as.numeric(CalledCount),y=derivedAlleles_rescaled))+
  geom_point()+
  facet_wrap(~category,scales="free")
########### plot just 0.8 ###########
# p5 <- ggplot(allCountsSynMis[allCountsSynMis$pop2!="BAJ" & allCountsSynMis$missingLevel==0.8,],aes(x=pop2,y=derivedAlleles_rescaled,fill=pop2))+
#   geom_boxplot()+
#   geom_point(alpha=0.4)+
#   facet_wrap(~category,scales="free")+
#   theme_bw()+
#   scale_fill_manual(values=c(colors['CA'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
#   theme(legend.position = "None")+
#   ggtitle("derived alleles (normalized)")+
#   theme(text=element_text(size=14))+
#   xlab("")
# p5
############ try counts with removing all fixed sites - maybe they are biasing things strangely #######
# this is with 0.8 missingness filter
newSynCountsExpt = read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/compareMissense_Syn/20181119/countsOfGenotypesPerIndividual/noFixedSites_minCallRate_0.8/syn.countsPerIndividual.minCall.0.8.countOfHomAltRefHet.txt",header=T)
head(newSynCountsExpt)
# merge with total callable
newSynCountsExpt_plusTotalCall <- merge(newSynCountsExpt,totalCallableCDSPerInd,by.x="individual",by.y="id")
head(newSynCountsExpt_plusTotalCall)
newSynCountsExpt_plusTotalCall$derivedAlleles <- 2*(newSynCountsExpt_plusTotalCall$HomAltCount)+newSynCountsExpt_plusTotalCall$HetCount
newSynCountsExpt_plusTotalCall$derivedAlleles_frac <- newSynCountsExpt_plusTotalCall$derivedAlleles/newSynCountsExpt_plusTotalCall$CalledCount
averageCalled = mean(totalCallableCDSPerInd$CalledCount)
newSynCountsExpt_plusTotalCall$derivedAlleles_rescaled <- newSynCountsExpt_plusTotalCall$derivedAlleles_frac * averageCalled

ggplot(newSynCountsExpt_plusTotalCall,aes(x=individual,y=derivedAlleles_rescaled))+
  geom_point()
ggplot(newSynCountsExpt_plusTotalCall,aes(x=CalledCount,y=derivedAlleles_rescaled))+
  geom_point()
