############ need to update this script to work on different missingness filters #######

require(ggplot2)
genotypeDate="20181119"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/compareMissense_Syn/",genotypeDate,"/countsOfGenotypesPerIndividual/",sep="")
colorPal=RColorBrewer::brewer.pal(n=8,name = "Dark2")
# make Baja brown
colors=list(CA=colorPal[1],BAJ=colorPal[7],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors
pops=c("CA","AK","AL","COM","KUR") # excluding BAJA because too low sample size
##### gather data ####
categories=c("syn","missense")
allCountsSynMis <- data.frame()
for(cat in categories){
counts <- read.table(paste(data.dir,cat,".countsPerIndividual.countOfHomAltRefHet.txt",sep=""),header=T) ### DONT USE! this has no missingness! use missingness filters (see below for example)
counts$derivedAlleles <- (2*counts$HomAltCount) + counts$HetCount # note -- Hom ref counts are very low because these are all snps, many of which are fixed 1/1 against ferret (elevates count of hom alt a lot); monomorphic ref sites will contain all the 0/0s.
counts$category <- cat
allCountsSynMis = rbind(allCountsSynMis,counts)
}
##### need to normalize by the total callable cds sites per individual #####
totalCallableCDSPerInd <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/compareMissense_Syn/20181119/countsOfGenotypesPerIndividual/callableCDSSitesPerIndividual.20191010.txt",header=T)

allCountsSynMis_plusCallable <- merge(allCountsSynMis,totalCallableCDSPerInd,by.x = "individual",by.y="id")
# get average called count across all inds:
avgCallableCDS = mean(totalCallableCDSPerInd$CalledCount)
# this should give those counts
# then what you want to do is:
allCountsSynMis_plusCallable$derivedAlleles_frac <- (allCountsSynMis_plusCallable$derivedAlleles / allCountsSynMis_plusCallable$CalledCount)
allCountsSynMis_plusCallable$derivedAlleles_normalized <- (allCountsSynMis_plusCallable$derivedAlleles / allCountsSynMis_plusCallable$CalledCount) * avgCallableCDS # normalize by dividing by total callable cds sites and mult by avg called sites (don't need to worry about sites vs alleles because factor of two factors out with site counts ratio)
allCountsSynMis_plusCallable$homAlt_normalized <- (allCountsSynMis_plusCallable$HomAltCount / allCountsSynMis_plusCallable$CalledCount) * avgCallableCDS # normalize by dividing by total callable cds sites and mult by avg called sites (don't need to worry about sites vs alleles because factor of two factors out with site counts ratio)

##### get populations: ####
allCountsSynMis_plusCallable$pop1 <- unlist(lapply(strsplit(as.character(allCountsSynMis_plusCallable$individual),"_"),"[",3))
# make Commanders one population :
allCountsSynMis_plusCallable$pop2 <- allCountsSynMis_plusCallable$pop1
allCountsSynMis_plusCallable[allCountsSynMis_plusCallable$pop1 %in% c("BER","MED"),]$pop2 <- "COM"

######## order factors ######

allCountsSynMis_plusCallable$category <- factor(allCountsSynMis_plusCallable$category,levels=c("syn","missense"))
allCountsSynMis_plusCallable$pop2 <- factor(allCountsSynMis_plusCallable$pop2,levels=c("CA","BAJ","AK","AL","COM","KUR"))
##### plot derived alleles (normalized): ####
p1a <- ggplot(allCountsSynMis_plusCallable,aes(x=pop2,y=derivedAlleles_normalized,fill=pop2))+
  geom_boxplot()+
  facet_wrap(~category,scales="free")+
  theme_bw()+
  scale_fill_manual(values=c(colors['CA'],colors['BAJ'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  ggtitle("Derived Alleles (normalized)")+
  theme(legend.position = "None")
p1a  
ggsave(paste(data.dir,"/derivedAlleles.Perpop.noMissingnessFilter.withBaja.pdf",sep=""),p1a,height=5,width=7)

#### drop Baja: ##### 
p1b <- ggplot(allCountsSynMis_plusCallable[allCountsSynMis_plusCallable$pop2!="BAJ",],aes(x=pop2,y=derivedAlleles_normalized,fill=pop2))+
  geom_boxplot()+
  geom_point(alpha=0.4)+
  facet_wrap(~category,scales="free")+
  theme_bw()+
  scale_fill_manual(values=c(colors['CA'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+   theme(legend.position = "none")+
  ggtitle("Derived Alleles (normalized)")
p1b
ggsave(paste(data.dir,"/derivedAlleles.Perpop.noMissingnessFilter.noBaja.pdf",sep=""),p1b,height=5,width=7)
###### plot homozygous alt #######
p2 <- ggplot(allCountsSynMis_plusCallable[allCountsSynMis_plusCallable$pop2!="BAJ",],aes(x=pop2,y=homAlt_normalized,fill=pop2))+
  geom_boxplot()+
  geom_point(alpha=0.4)+
  facet_wrap(~category,scales="free")+
  theme_bw()+
  scale_fill_manual(values=c(colors['CA'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme(legend.position = "None")+
  ggtitle("Homozygous derived genotypes (normalized)")
p2
ggsave(paste(data.dir,"/homAltGTs.Perpop.noMissingnessFilter.noBaja.pdf",sep=""),p2,height=5,width=7)

### idea: recalculate removing sites with excess missing data? ###
### removed all sites with any missing data. so then I think don't have to rescale because all sites would have the same number of cds sites with no missing data (would have to rescale if I go down to 80%)
# see how it's looking
# this is syn. only:
syntest <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/troubleshooting/sandboxCountDerived/noMissingDataAllowed/syn.countsPerIndividual.NoMissing.countOfHomAltRefHet.txt",header=T)
##### get populations: ####
syntest$pop1 <- unlist(lapply(strsplit(as.character(syntest$individual),"_"),"[",3))
# make Commanders one population :
syntest$pop2 <- syntest$pop1
syntest[syntest$pop1 %in% c("BER","MED"),]$pop2 <- "COM"
head(syntest)
syntest$derivedAlleles <- (2*syntest$HomAltCount) + syntest$HetCount
ggplot(syntest,aes(x=pop2,y=derivedAlleles,fill=pop2))+
  geom_boxplot()+
  scale_fill_manual(values=c(colors['CA'],colors['BAJ'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme_bw()
  
ggplot(syntest,aes(x=pop2,y=HomAltCount,fill=pop2))+
  geom_boxplot()+
  scale_fill_manual(values=c(colors['CA'],colors['BAJ'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme_bw()

ggplot(syntest,aes(x=pop2,y=HetCount,fill=pop2))+
  geom_boxplot()+
  scale_fill_manual(values=c(colors['CA'],colors['BAJ'],colors['AK'],colors['AL'],colors['COM'],colors['KUR']))+
  theme_bw()
############## DO SOME EXPLORING ##########
ggplot(allCountsSynMis_plusCallable,aes(x=CalledCount,y=derivedAlleles_frac,color=pop2))+
  geom_point()+
  facet_wrap(~category,scales="free")+
  theme_bw()+
  ggtitle("Fraction of derived alleles is correlated\nwith number of called cds sites (drop out isn't random?)")
ggplot(allCountsSynMis_plusCallable,aes(x=CalledCount,y=homAlt_normalized,color=pop2))+
  geom_point()+
  facet_wrap(~category,scales="free")+
  theme_bw()+
  ggtitle("Number of hom-alt sites is correlated\nwith number of called cds sites (drop out isn't random?)")

## look at this with missingness filter of 0.8 -- might be a lot better.

################## explore neutral regions (need to normalize by neutral _all) #####
neutral <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/compareMissense_Syn/20181119/countsOfGenotypesPerIndividual/neutral.countsPerIndividual.countOfHomAltRefHet.txt",header=T)
head(neutral)
neutral$derivedAlleles <- (2*neutral$HomAltCount) + neutral$HetCount
ggplot(neutral,aes(x=individual,y=derivedAlleles))+
  geom_point()
# normalize?

# maybe use neutral to get a neutral drop out rate to use in fudge factor! <- cool idea
# normalize by total called neut sites first though