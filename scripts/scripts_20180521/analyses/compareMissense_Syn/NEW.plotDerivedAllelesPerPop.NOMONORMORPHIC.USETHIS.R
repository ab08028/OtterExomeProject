require(ggplot2)
genotypeDate="20181119"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/compareMissense_Syn/",genotypeDate,"/countsOfGenotypesPerIndividual/",sep="")
colorPal=RColorBrewer::brewer.pal(n=8,name = "Dark2")
# make Baja brown
colors=list(CA=colorPal[1],BAJ=colorPal[7],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5]) # your population colors
pops=c("CA","AK","AL","COM","KUR") # excluding BAJA because too low sample size
date="20191014" # annoying -- maybe don't output this in future
categories=c("syn","missense")
missingness=c("0.90","0.95") # this is max missingness threshold from vcftools which is unhelpfully named. it's really the min call rate, so 0 refers to any amount of missingness being allowed  (min call rate is 0); 0.8 allows 20% missingness; 1 allows no missingness. very confusing nomenclature, I know. Leaving 0 out because already plotted using old script and it has a slightly different filenames structure -- look in "olderScripts" to see plottin g of missingness  = 0 (all missing allowed)
allCountsSynMis <- data.frame()
for(missing in missingness){
  for(cat in categories){
    counts <- read.table(paste(data.dir,"/minCallRate_",missing,".NOMONOMORPHIC/",cat,".countsPerIndividual.minCall.",missing,".NOMONOMORPHIC.countOfHomAltRefHet.txt",sep=""),header=T) ### DONT USE! this has no missingness! use missingness filters (see below for example)
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
ggsave(paste(data.dir,"/derivedAlleles.Perpop.MultMissingFilts.withBaja.NOMONOMORPHIC.pdf",sep=""),p1,height=7,width=8)
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
ggsave(paste(data.dir,"/derivedAlleles.Perpop.MultMissingFilts.noBaja.NOMONOMORPHIC.pdf",sep=""),p2,height=7,width=8)

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
ggsave(paste(data.dir,"/homAlt.Perpop.MultMissingFilts.noBaja.NOMONOMORPHIC.pdf",sep=""),p3,height=7,width=8)
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
ggsave(paste(data.dir,"/homAlt.Perpop.MultMissingFilts.NOMONOMORPHIC.pdf",sep=""),p4,height=15,width=8)

#### plot hets, no baja: NOTE: lack of monomorphic sites doesn't make a diff for hets
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
ggsave(paste(data.dir,"/hets.Perpop.MultMissingFilts.noBaja.NOMONOMORPHIC.pdf",sep=""),p4b,height=7,width=8)
