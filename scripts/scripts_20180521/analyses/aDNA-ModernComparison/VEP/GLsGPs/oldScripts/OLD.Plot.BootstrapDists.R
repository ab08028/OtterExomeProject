require(ggplot2)
require(reshape2)
require(dplyr)
###### Plot distribution
#dates=c("20190701-lowcov-AFprior-MajorMinor4", "20190701-highcov-AFprior-MajorMinor4")

# need ind file for HC LC to identify individuals as ancient/modern
# need point estimates to put in as points
#  
minDepth=2
minGP=0.95
minInd=1

###################### LOW COVERAGE ############
date="20190701-lowcov-AFprior-MajorMinor4"
bamListDir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/"
# these point estimates were old and are not based on windowing the genome; if all windows are included should be the same, but I want to use the point estimate based on my good windows and the averaging of windows.
indir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/compareMisSynDists_withBootstraps/",date,"/PointEstsPlusBootstraps/",sep="")
#pointDir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/sumGPsGLsPerVEPCategory/",date,"/",sep="")
#bootsDir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/compareMisSynDists_withBootstraps/",date,"/",sep="")
bamListFile="SampleIDsInOrder.LowCoverageOnly.BeCarefulOfOrder.txt" ### CHECK CAREFULLY
bamList=read.table(paste(bamListDir,bamListFile,sep=""),header=F)
plot.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/plots/",date,"/",sep="")
dir.create(plot.dir,showWarnings = F)
# this is the order of samples in low coverage (is in diff order for high coverage)
# so ind0 = bamList[1,]
# point estimates:
# point estimates are gotten from summing up GPs across all rescaled genome windows that pass filters (within each window, ancient and modern were rescaled by avg called cds sites in that window)

pointEsts <- read.table(paste(indir,"angsdOut.mappedTomfur.Bootstraps.GPs.ProbCutoff.",minGP,".DepthCutoff.",minDepth,".minInd.",minInd,".",date,".Modern.Ancient.PointEstimatesBasedonGoodBins.txt",sep=""),header = T)
# then want to plot these point estimates on top of distributions as points (don't freak out if it's weird -- this is all new and you're trouble shooting)

allIndBootstraps=data.frame()
for(Ind in seq(0,8)){
  input <- read.table(paste(bootsDir,"/angsdOut.mappedTomfur.Bootstraps.GPs.ProbCutoff.",minGP,".DepthCutoff.",minDepth,".minInd.",minInd,".",date,".Ind.",Ind,".allBoots.txt",sep=""),header=T) # eventually change this to reflect filters
  input$individualID <- bamList[Ind+1,] # need to add 1 because bamList indices are 1 based and Ind numbers are 0 based; so Ind0 is in position 1 of bamList, Ind2 is in position 3, etc.
  input$minDepth <- minDepth
  input$minGP <- minGP
  input$minInd <- minInd
  allIndBootstraps <- rbind(allIndBootstraps,input)
}

######## assign ancient/modern designations #######
allIndBootstraps$ancOrModern <- "Modern"
allIndBootstraps[grepl("^A",allIndBootstraps$individualID),]$ancOrModern <- "Ancient"

########  get derived alleles : ###########
allIndBootstraps$DerivedAlleles <- (2*allIndBootstraps$homAlt)+allIndBootstraps$het
allIndBootstraps$DerivedAllelesFrac <- allIndBootstraps$DerivedAlleles / allIndBootstraps$totalCDSSitesInBoot
meanCallableSites = unique(pointEsts$meanCallableSitesAllSitesAcrossInds)
allIndBootstraps$DerivedAlleles_Rescaled <- allIndBootstraps$DerivedAllelesFrac * meanCallableSites


########## reorder factors: #########
allIndBootstraps$category <- factor(allIndBootstraps$category,levels=c("synonymous","missense","stop_gained"))

######## get point estimates to be dots inside dists ############
test = pointEsts[,c("sample","callableSites_AllCDSSites","sumHetGLsOrGPs","sumHetGLsOrGPs_TransversionsOnly","sumHomAltGLsOrGPs","sumHomAltGLsOrGPs_TransversionsOnly", "sumHomRefGLsOrGPs","category", "meanCallableSitesAllSitesAcrossInds","ancOrModern","derivedAllelesAll","derivedAllelesAll_TV","derivedAllelesAll_meanScaled","derivedAllelesAll_TV_meanScaled","HomAlt_meanScaled","HomAlt_meanScaled_TV")]
test2=pointEsts[,c("sample","category", "ancOrModern","derivedAllelesAll_meanScaled","derivedAllelesAll_TV_meanScaled")]
test2melt <- melt(test2)
# label as Tv only:
test2melt$siteType <- "Ti+Tv"
test2melt[grepl("TV",test2melt$variable),]$siteType <- "TvOnly"
test2melt$category <- factor(test2melt$category,levels=c("synonymous","missense","stop_gained"))

############ plot derived alleles per individual ###############
p1a <- ggplot(allIndBootstraps,aes(x=individualID,y=DerivedAlleles_Rescaled))+
  geom_violin()+
  geom_point(data=test2melt[test2melt$sample %in% allIndBootstraps$individualID,],aes(x=sample,y=value))+
  facet_wrap(siteType~category,scales="free")+
  theme_bw()+
  ggtitle(paste("**Dowsampled Modern**\nDerived Alleles per individual\n",date,"\nMinGP:",minGP,";minDepth:",minDepth,"minInd:",minInd,sep=""))

p1a

ggsave(paste(plot.dir,"derivedAlleles.PerIndividual.minGP.",minGP,".minDepth.",minDepth,".minInd",minInd,".pdf",sep=""),p1a,device="pdf",height=5,width=7)

############ plot derived alleles per anc/modern ###############
p1b <- ggplot(allIndBootstraps,aes(x=ancOrModern,y=DerivedAlleles_Rescaled))+
  geom_violin(aes(fill=ancOrModern),alpha=0.8)+
  geom_point(data=test2melt[test2melt$sample %in% allIndBootstraps$individualID,],aes(x=ancOrModern,y=value))+
  facet_wrap(siteType~category,scales="free")+
  theme_bw()+
  scale_fill_manual(values=c("dodgerblue","orange"))+
  ggtitle(paste("**Dowsampled Modern**\nDerived Alleles compared between ancient and modern\n",date,"\nMinGP:",minGP,";minDepth:",minDepth,"minInd:",minInd,sep=""))

p1b

ggsave(paste(plot.dir,"derivedAlleles.Anc.Modern.minGP.",minGP,".minDepth.",minDepth,".minInd",minInd,".pdf",sep=""),p1b,device="pdf",height=5,width=7)

############### plot homozygous alternate #############


# plot hom derived:
# just the homAlt stat:homAlt
test3=pointEsts[,c("sample","category", "ancOrModern","HomAlt_meanScaled","HomAlt_meanScaled_TV")]
test3melt <- melt(test3)
# label as Tv only:
test3melt$siteType <- "Ti+Tv"
test3melt[grepl("TV",test3melt$variable),]$siteType <- "TvOnly"
test3melt$category <- factor(test3melt$category,levels=c("synonymous","missense","stop_gained"))


p2a <- ggplot(allIndBootstraps,aes(x=individualID,y=homAlt))+
  geom_violin()+
  geom_point(data=test3melt[test3melt$sample %in% allIndBootstraps$individualID,],aes(x=sample,y=value))+
  facet_wrap(siteType~category,scales="free")+
  theme_bw()+
  ggtitle(paste("**Dowsampled Modern**\nHomozygous alternate genotypes per individual/modern\n",date,"\nMinGP:",minGP,";minDepth:",minDepth,"minInd:",minInd,sep=""))
  p2a

ggsave(paste(plot.dir,"homAltGenotypes.PerIndividual.minGP.",minGP,".minDepth.",minDepth,".minInd",minInd,".pdf",sep=""),p2a,device="pdf",height=5,width=7)


p2b <- ggplot(allIndBootstraps,aes(x=ancOrModern,y=homAlt))+
  geom_violin(aes(fill=ancOrModern),alpha=0.8)+
  geom_point(data=test3melt[test3melt$sample %in% allIndBootstraps$individualID,],aes(x=ancOrModern,y=value))+
  facet_wrap(siteType~category,scales="free")+
  theme_bw()+
  scale_fill_manual(values=c("dodgerblue","orange"))+
  ggtitle(paste("**Dowsampled Modern**\nHomozygous alternate genotypes compared between ancient/modern\n",date,"\nMinGP:",minGP,";minDepth:",minDepth,"minInd:",minInd,sep=""))+
  ylab("Homozygous alternate genotypes (rescaled)")+
  xlab("")
p2b

ggsave(paste(plot.dir,"homAltGenotypes.Anc.Modern.minGP.",minGP,".minDepth.",minDepth,".minInd",minInd,".pdf",sep=""),p2b,device="pdf",height=5,width=7)

##### trouble shooting ############
# fix stopgained vs stop_gained


# look at pt estimates from bins: for Ind.2 only:
sumsPerBin <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/sumGPsGLsPerVEPCategory/20190701-lowcov-AFprior-MajorMinor4/bootstraps/angsdOut.mappedTomfur.superfile.GPs.Ind.2.sumsPerBin.txt",header=T)
require(dplyr)
missenseOnly = sumsPerBin[grep("missense_variant",sumsPerBin$Consequence),]

missenseOnly %>%
  group_by(sites) %>%
  summarise(totalHomRef= sum(sumHomRef),totalHomAlt=sum(sumHomAlt),totalHet=sum(sumHet))

#View(pointEsts)
#View(missenseOnly)




###################### HIGH COVERAGE ############
date="20190701-highcov-AFprior-MajorMinor4"
bamListDir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_calling_aDNA/bamLists/"
pointDir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/sumGPsGLsPerVEPCategory/",date,"/",sep="")
bootsDir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/compareMisSynDists_withBootstraps/",date,"/",sep="")
bamListFile="SampleIDsInOrder.HighCoverageAndADNAOnly.BeCarefulOfOrder.txt" ### CHECK CAREFULLY
bamList=read.table(paste(bamListDir,bamListFile,sep=""),header=F)
plot.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/plots/",date,"/",sep="")
dir.create(plot.dir,showWarnings = F)
# this is the order of samples in low coverage (is in diff order for high coverage)
# so ind0 = bamList[1,]
# point estimates:
# point estimates gotten from summing up GPs acrosss sites passing filters and rescaling by teh avg called cds sites across all individuals (done in Step 5)
pointEsts <- read.table(paste(pointDir,"angsdOut.mappedTomfur.hetHomTotals.GPs.ProbCutoff.",minGP,".DepthCutoff.",minDepth,".minInd.",minInd,".",date,".AllPointEstimates.allCategories.txt",sep=""),header = T)
# make stopgained stop_gained to match
pointEsts$category2 <- as.character(pointEsts$category)
pointEsts[pointEsts$category=="stopgained",]$category2 <- "stop_gained"
pointEsts$category <- pointEsts$category2
# then want to plot these point estimates on top of distributions as points (don't freak out if it's weird -- this is all new and you're trouble shooting)

allIndBootstraps=data.frame()
for(Ind in seq(0,8)){
  input <- read.table(paste(bootsDir,"/angsdOut.mappedTomfur.Bootstraps.GPs.ProbCutoff.",minGP,".DepthCutoff.",minDepth,".minInd.",minInd,".",date,".Ind.",Ind,".allBoots.txt",sep=""),header=T) # eventually change this to reflect filters
  input$individualID <- bamList[Ind+1,] # need to add 1 because bamList indices are 1 based and Ind numbers are 0 based; so Ind0 is in position 1 of bamList, Ind2 is in position 3, etc.
  input$minDepth <- minDepth
  input$minGP <- minGP
  input$minInd <- minInd
  allIndBootstraps <- rbind(allIndBootstraps,input)
}

######## assign ancient/modern designations #######
allIndBootstraps$ancOrModern <- "Modern"
allIndBootstraps[grepl("^A",allIndBootstraps$individualID),]$ancOrModern <- "Ancient"

########  get derived alleles : ###########
allIndBootstraps$DerivedAlleles <- (2*allIndBootstraps$homAlt)+allIndBootstraps$het
allIndBootstraps$DerivedAllelesFrac <- allIndBootstraps$DerivedAlleles / allIndBootstraps$totalCDSSitesInBoot
meanCallableSites = unique(pointEsts$meanCallableSitesAllSitesAcrossInds)
allIndBootstraps$DerivedAlleles_Rescaled <- allIndBootstraps$DerivedAllelesFrac * meanCallableSites


########## reorder factors: #########
allIndBootstraps$category <- factor(allIndBootstraps$category,levels=c("synonymous","missense","stop_gained"))

######## get point estimates to be dots inside dists ############
test = pointEsts[,c("sample","callableSites_AllCDSSites","sumHetGLsOrGPs","sumHetGLsOrGPs_TransversionsOnly","sumHomAltGLsOrGPs","sumHomAltGLsOrGPs_TransversionsOnly", "sumHomRefGLsOrGPs","category", "meanCallableSitesAllSitesAcrossInds","ancOrModern","derivedAllelesAll","derivedAllelesAll_TV","derivedAllelesAll_meanScaled","derivedAllelesAll_TV_meanScaled","HomAlt_meanScaled","HomAlt_meanScaled_TV")]
test2=pointEsts[,c("sample","category", "ancOrModern","derivedAllelesAll_meanScaled","derivedAllelesAll_TV_meanScaled")]
test2melt <- melt(test2)
# label as Tv only:
test2melt$siteType <- "Ti+Tv"
test2melt[grepl("TV",test2melt$variable),]$siteType <- "TvOnly"
test2melt$category <- factor(test2melt$category,levels=c("synonymous","missense","stop_gained"))
############ plot derived alleles per individual ###############
p3a <- ggplot(allIndBootstraps,aes(x=individualID,y=DerivedAlleles_Rescaled))+
  geom_violin()+
  geom_point(data=test2melt[test2melt$sample %in% allIndBootstraps$individualID,],aes(x=sample,y=value))+
  facet_wrap(siteType~category,scales="free")+
  theme_bw()+
  ggtitle(paste("**High coverage Modern**\nDerived Alleles per individual\n",date,"\nMinGP:",minGP,";minDepth:",minDepth,"minInd:",minInd,sep=""))

p3a

ggsave(paste(plot.dir,"derivedAlleles.PerIndividual.minGP.",minGP,".minDepth.",minDepth,".minInd",minInd,".pdf",sep=""),p3a,device="pdf",height=5,width=7)

############ plot derived alleles per anc/modern ###############
p3b <- ggplot(allIndBootstraps,aes(x=ancOrModern,y=DerivedAlleles_Rescaled))+
  geom_violin(aes(fill=ancOrModern),alpha=0.8)+
  geom_point(data=test2melt[test2melt$sample %in% allIndBootstraps$individualID,],aes(x=ancOrModern,y=value))+
  facet_wrap(siteType~category,scales="free")+
  theme_bw()+
  scale_fill_manual(values=c("dodgerblue","orange"))+
  ggtitle(paste("**High coverage Modern**\nDerived Alleles compared between ancient and modern\n",date,"\nMinGP:",minGP,";minDepth:",minDepth,"minInd:",minInd,sep=""))

p3b

ggsave(paste(plot.dir,"derivedAlleles.Anc.Modern.minGP.",minGP,".minDepth.",minDepth,".minInd",minInd,".pdf",sep=""),p3b,device="pdf",height=5,width=7)

############### plot homozygous alternate #############


# plot hom derived:
# just the homAlt stat:homAlt
test3=pointEsts[,c("sample","category", "ancOrModern","HomAlt_meanScaled","HomAlt_meanScaled_TV")]
test3melt <- melt(test3)
# label as Tv only:
test3melt$siteType <- "Ti+Tv"
test3melt[grepl("TV",test3melt$variable),]$siteType <- "TvOnly"


p4a <- ggplot(allIndBootstraps,aes(x=individualID,y=homAlt))+
  geom_violin()+
  geom_point(data=test3melt[test3melt$sample %in% allIndBootstraps$individualID,],aes(x=sample,y=value))+
  facet_wrap(siteType~category,scales="free")+
  theme_bw()+
  ggtitle(paste("**High coverage Modern**\nHomozygous alternate genotypes per individual/modern\n",date,"\nMinGP:",minGP,";minDepth:",minDepth,"minInd:",minInd,sep=""))
p4a

ggsave(paste(plot.dir,"homAltGenotypes.PerIndividual.minGP.",minGP,".minDepth.",minDepth,".minInd",minInd,".pdf",sep=""),p4a,device="pdf",height=5,width=7)


p4b <- ggplot(allIndBootstraps,aes(x=ancOrModern,y=homAlt))+
  geom_violin(aes(fill=ancOrModern),alpha=0.8)+
  geom_point(data=test3melt[test3melt$sample %in% allIndBootstraps$individualID,],aes(x=ancOrModern,y=value))+
  facet_wrap(siteType~category,scales="free")+
  theme_bw()+
  scale_fill_manual(values=c("dodgerblue","orange"))+
  ggtitle(paste("**High coverage Modern**\nHomozygous alternate genotypes compared between ancient/modern\n",date,"\nMinGP:",minGP,";minDepth:",minDepth,"minInd:",minInd,sep=""))+
  ylab("Homozygous alternate genotypes (rescaled)")+
  xlab("")
p4b

ggsave(paste(plot.dir,"homAltGenotypes.Anc.Modern.minGP.",minGP,".minDepth.",minDepth,".minInd",minInd,".pdf",sep=""),p4b,device="pdf",height=5,width=7)

