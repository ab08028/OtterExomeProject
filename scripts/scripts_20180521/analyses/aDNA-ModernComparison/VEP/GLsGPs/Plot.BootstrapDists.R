require(ggplot2)
require(reshape2)
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
pointDir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/sumGPsGLsPerVEPCategory/",date,"/",sep="")
bootstrapDir=paste(pointDir,"bootstraps/",sep="")
bamListFile="SampleIDsInOrder.LowCoverageOnly.BeCarefulOfOrder.txt" ### CHECK CAREFULLY
bamList=read.table(paste(bamListDir,bamListFile,sep=""),header=F)
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
for(Ind in c(2,4,6)){
  input <- read.table(paste(bootstrapDir,"/angsdOut.mappedTomfur.superfile.GPs.Ind.",Ind,".allBoots.txt",sep=""),header=T) # eventually change this to reflect filters
  input$individualID <- bamList[Ind+1,] # need to add 1 because bamList indices are 1 based and Ind numbers are 0 based; so Ind0 is in position 1 of bamList, Ind2 is in position 3, etc.
  input$minDepth <- minDepth
  input$minGP <- minGP
  input$minInd <- minInd
  allIndBootstraps <- rbind(allIndBootstraps,input)
}

# get derived alleles:
allIndBootstraps$DerivedAlleles <- (2*allIndBootstraps$homAlt)+allIndBootstraps$het
allIndBootstraps$DerivedAllelesFrac <- allIndBootstraps$DerivedAlleles / allIndBootstraps$totalCDSSitesInBoot
meanCallableSites = unique(pointEsts$meanCallableSitesAllSitesAcrossInds)
allIndBootstraps$DerivedAlleles_Rescaled <- allIndBootstraps$DerivedAllelesFrac * meanCallableSites

# bootstraps
test = pointEsts[,c("sample","callableSites_AllCDSSites","sumHetGLsOrGPs","sumHetGLsOrGPs_TransversionsOnly","sumHomAltGLsOrGPs","sumHomAltGLsOrGPs_TransversionsOnly", "sumHomRefGLsOrGPs","category", "meanCallableSitesAllSitesAcrossInds","ancOrModern","derivedAllelesAll","derivedAllelesAll_TV","derivedAllelesAll_meanScaled","derivedAllelesAll_TV_meanScaled","HomAlt_meanScaled","HomAlt_meanScaled_TV")]
test2=pointEsts[,c("sample","category", "ancOrModern","derivedAllelesAll_meanScaled","derivedAllelesAll_TV_meanScaled")]
test2melt <- melt(test2)
# label as Tv only:
test2melt$siteType <- "Ti+Tv"
test2melt[grepl("TV",test2melt$variable),]$siteType <- "TvOnly"
ggplot(allIndBootstraps,aes(x=individualID,y=DerivedAlleles_Rescaled))+
  geom_violin()+
  geom_point(data=test2melt[test2melt$sample %in% allIndBootstraps$individualID,],aes(x=sample,y=value))+
  facet_wrap(siteType~category,scales="free")

# test:
allIndBootstraps$DerivedAlleles <- (2*allIndBootstraps$homAlt)+allIndBootstraps$het

##### weird things to fix:
# fix stopgained vs stop_gained


# look at pt estimates from bins: for Ind.2 only:
sumsPerBin <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/sumGPsGLsPerVEPCategory/20190701-lowcov-AFprior-MajorMinor4/bootstraps/angsdOut.mappedTomfur.superfile.GPs.Ind.2.sumsPerBin.txt",header=T)
require(dplyr)
missenseOnly = sumsPerBin[grep("missense_variant",sumsPerBin$Consequence),]

missenseOnly %>%
  group_by(sites) %>%
  summarise(totalHomRef= sum(sumHomRef),totalHomAlt=sum(sumHomAlt),totalHet=sum(sumHet))

View(pointEsts)
View(missenseOnly)
