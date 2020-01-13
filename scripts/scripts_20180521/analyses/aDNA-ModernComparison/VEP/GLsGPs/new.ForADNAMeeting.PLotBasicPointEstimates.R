############## Go back before things were rescaled per window
# just do something way more simple ala orlando (who didn't do any bootstrapping)
# and show hom-alt GTs in ancient and modern rescaled by total hom-alt genotypes

### so going toward way less transformed data
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/sumGPsGLsPerVEPCategory_notBasedOnWindows/"
plot.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/plots/basicPointEstimates_noWindowsOrWIndowFiltering_foraDNAPresentation/"
dates=c("20190701-highcov-AFprior-MajorMinor4","20190701-lowcov-AFprior-MajorMinor4")
depth=2
categories=c("synonymous","missense","stopgained","CDS")
allInputs = data.frame()
for(angsdDate in dates){
  # get total cds sites:
  cat="CDS"
  CDSinput <- read.table(paste(data.dir,angsdDate,"/angsdOut.mappedTomfur.hetHomTotals.GPs.ProbCutoff.0.95.DepthCutoff.",depth,".minInd.1.",angsdDate,".",cat,".txt",sep=""),header=T)
  head(CDSinput)
  totals <- CDSinput %>%
    filter(Filter_PerIndividualDepthMinimum==depth) %>% # this is just a back up check
    group_by(sample) %>%
    summarise(totalCallableCDS=sum(callableSites),totalCDSHomGTS=sum(sum(sumHomAltGLsOrGPs),sum(sumHomRefGLsOrGPs)))
  meanCallableCDS <- mean(totals$totalCallableCDS)
  meanHomozygousCDS <- mean(totals$totalCDSHomGTS)
  for(cat in categories){
    ####### want to plot the point estimate (don't worry about boots right now)
    # of syn mis and sg sites per individual
    # rescaled either by totall called cds sites, or by cds hom-gts (what orlando does to apparently control for inbreeding somehow)
    # /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/sumGPsGLsPerVEPCategory_notBasedOnWindows/20190701-lowcov-AFprior-MajorMinor4/angsdOut.mappedTomfur.hetHomTotals.GPs.ProbCutoff.0.95.DepthCutoff.2.minInd.1.20190701-lowcov-AFprior-MajorMinor4.synonymous.txt
    
    input <- read.table(paste(data.dir,angsdDate,"/angsdOut.mappedTomfur.hetHomTotals.GPs.ProbCutoff.0.95.DepthCutoff.",depth,".minInd.1.",angsdDate,".",cat,".txt",sep=""),header=T)
    inputPlustTotals <- merge(input,totals,by="sample")
    inputPlustTotals$frac_HomAltGLsOrGPs_TVOnly_FracScaledByTotalHomGTs <- inputPlustTotals$sumHomAltGLsOrGPs_TransversionsOnly / inputPlustTotals$totalCDSHomGTS
    inputPlustTotals$category <- cat
    inputPlustTotals$date <- angsdDate
    inputPlustTotals$meanCallableCDS_allInds <- meanCallableCDS
    inputPlustTotals$meanCallableCDS_allInds_Homozygotes <- meanHomozygousCDS
    allInputs <- rbind(allInputs,inputPlustTotals)
    
  }}
### just hc 
allInputs$group <- "modern"
allInputs[grep("^A",allInputs$sample),]$group <- "ancient"
# exclude cds
#order factors:
allInputs$category <- factor(allInputs$category, levels=c("synonymous","missense","stopgained","CDS"))

## want a common multiplier for the hc and lc that is the same -- let's just do 10M sites to get a sense (so we can compare directly)

p1 <- ggplot(allInputs[allInputs$date=="20190701-highcov-AFprior-MajorMinor4" & allInputs$category!="CDS",],aes(x=group,y=frac_HomAltGLsOrGPs_TVOnly_FracScaledByTotalHomGTs))+
  geom_boxplot(alpha=0.8,aes(fill=group))+
  geom_point()+
  #geom_jitter()+
  facet_wrap(~category,scales="free")+
  ylab("Fraction of homozygous genotypes that are hom-alt TVs")+
  theme_bw()+
  ggtitle("Fraction of homozygous genotypes that are hom-alt transversions\nNon-Downsampled Modern Data")+
  theme(legend.position = "none")
p1
ggsave(paste(plot.dir,"fracHomGTsThatAreHomAlt.HighCov.pdf",sep=""),p1,height=5,width=7)


p2 <- ggplot(allInputs[allInputs$date=="20190701-lowcov-AFprior-MajorMinor4" & allInputs$category!="CDS",],aes(x=group,y=frac_HomAltGLsOrGPs_TVOnly_FracScaledByTotalHomGTs))+
  geom_boxplot(alpha=0.8,aes(fill=group))+
  geom_point()+
  #geom_jitter()+
  facet_wrap(~category,scales="free")+
  ylab("Fraction of homozygous genotypes that are hom-alt TVs")+
  theme_bw()+
  ggtitle("Fraction of homozygous genotypes that are hom-alt transversions\nDownsampled Modern Data")+
  theme(legend.position = "none")
p2
ggsave(paste(plot.dir,"fracHomGTsThatAreHomAlt.LowCov.pdf",sep=""),p2,height=5,width=7)


####### theoretically scale up fraction by 10M sites 
unifiedScalingFactor=10000000 # 10M sites
p3 <- ggplot(allInputs[allInputs$date=="20190701-lowcov-AFprior-MajorMinor4" & allInputs$category!="CDS",],aes(x=group,y=frac_HomAltGLsOrGPs_TVOnly_FracScaledByTotalHomGTs*totalCDSHomGTS))+
  geom_boxplot(alpha=0.8,aes(fill=group))+
  geom_point()+
  #geom_jitter()+
  facet_wrap(~category,scales="free")+
  ylab("Fraction of homozygous genotypes that are hom-alt TVs")+
  theme_bw()+
  ggtitle("Fraction of homozygous genotypes that are hom-alt transversions\nDownsampled Modern Data")+
  theme(legend.position = "none")
p3
ggsave(paste(plot.dir,"scaledUpBy10MHomozygousGenotypes.AsExample.NotRealScaling.pdf",sep=""),p3,height=5,width=7) 


######### just plot total number rescaled by average #######
p4 <- ggplot(allInputs[allInputs$date=="20190701-lowcov-AFprior-MajorMinor4" & allInputs$category!="CDS",],aes(x=group,y=frac_HomAltGLsOrGPs_TVOnly_FracScaledByTotalHomGTs*meanCallableCDS_allInds_Homozygotes))+
  geom_boxplot(alpha=0.8,aes(fill=group))+
  geom_point()+
  #geom_jitter()+
  facet_wrap(~category,scales="free")+
  ylab("homozygous-alt TV genotypes")+
  theme_bw()+
  ggtitle("homozygous-alt transversions genotypes\n(rescaled by avg called homozygous genotypes)\nDownsampled Modern Data")+
  theme(legend.position = "none")
p4
ggsave(paste(plot.dir,"scaledUpByMeanHomozygousGenotypes.pdf",sep=""),p4,height=5,width=7) 

allInputs$sampleSHORT <- unlist(lapply(strsplit(as.character(allInputs$sample),"_"),"[",1))
p5 <-  ggplot(allInputs[allInputs$date=="20190701-lowcov-AFprior-MajorMinor4" & allInputs$category!="CDS",],aes(x=sampleSHORT,y=frac_HomAltGLsOrGPs_TVOnly_FracScaledByTotalHomGTs*meanCallableCDS_allInds_Homozygotes,color=group))+
  geom_point()+
  #geom_jitter()+
  facet_wrap(~category,scales="free")+
  ylab("homozygous-alt TV genotypes")+
  theme_bw()+
  ggtitle("homozygous-alt transversions genotypes\n(rescaled by avg called homozygous genotypes)\nDownsampled Modern Data")+
  theme(legend.position = "none")+
  xlab("")
p5
