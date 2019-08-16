require(ggplot2)
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/sumGPsGLsPerVEPCategory/"

dates=c("20190701-highcov-AFprior-MajorMinor4","20190701-lowcov-AFprior-MajorMinor4")
minInds=c(1)
ProbCutoffs=c(0.95)
DepthCutoffs=c(1,2,4)
#DepthCutoffs=2
# cds is all coding; missense; synonymous (skipping SG for now since it's a bit weird after using --pick.)
categories=c("missense","synonymous","stopgained") # something's a bit odd with stopgained, too few? or too many in otter paper perhaps. Maybe --pick reclassifies them, that is most likely. then that's okay. because then it's a parameter choice and it is the same between each otter in genome paper, so even if absolute # is disputable, the *pattern* is the same. Okay. cool.

ref="mfur"
type="GPs"
for(date in dates){
  allInputsAllFilters <- data.frame()
  for(minInd in minInds){
    for(ProbCutoff in ProbCutoffs){
      for(DepthCutoff in DepthCutoffs){
        # get total cds sites:
        CDSinfile=paste(data.dir,date,"/angsdOut.mappedTo",ref,".hetHomTotals.",type,".ProbCutoff.",ProbCutoff,".DepthCutoff.",DepthCutoff,".minInd.",minInd,".",date,".CDS.txt",sep="")
        # out file name:
        #outfileName=paste(data.dir,date,"/angsdOut.mappedTo",ref,".hetHomTotals.",type,".ProbCutoff.",ProbCutoff,".DepthCutoff.",DepthCutoff,".minInd.",minInd,".",date,".AllPointEstimates.allCategories.txt",sep="")
        cdsSites <- read.table(CDSinfile,header=T,sep="\t")
        callableSites <- cdsSites[,c("sample","callableSites")]
        meanCallableSitesAcrossInds <- mean(callableSites$callableSites)
        # merge this with each input
        allCategories=data.frame()
        for(category in categories){
          infile=paste(data.dir,date,"/angsdOut.mappedTo",ref,".hetHomTotals.",type,".ProbCutoff.",ProbCutoff,".DepthCutoff.",DepthCutoff,".minInd.",minInd,".",date,".",category,".txt",sep="")
          input <- read.table(infile, header=T,sep="\t")
          input$category <- category
          #input$DepthCutoff <- DepthCutoff
          #input$ProbCutoff <- ProbCutoff
          #input$minInd <- minInd
          input$meanCallableSitesAllSitesAcrossInds <- meanCallableSitesAcrossInds
          input$date <- date
          input$ancOrModern <- "Modern"
          input[grep("^A",input$sample),]$ancOrModern <- "Ancient"
          comboInput <- merge(callableSites,input,by="sample",suffixes=c("_AllCDSSites"))
          allCategories <- rbind(allCategories,comboInput)
        }
        # rescale all categories by avg called sites
        allCategories$derivedAllelesAll <- allCategories$sumHetGLsOrGPs + 2*allCategories$sumHomAltGLsOrGPs
        allCategories$derivedAllelesAll_frac <- allCategories$derivedAllelesAll / allCategories$callableSites_AllCDSSites
        allCategories$derivedAllelesAll_meanScaled <- allCategories$derivedAllelesAll_frac * allCategories$meanCallableSitesAllSitesAcrossInds # don't need to do factor of 2 bc cancels out
        # transversions:
        allCategories$derivedAllelesAll_TV <- allCategories$sumHetGLsOrGPs_TransversionsOnly + 2*allCategories$sumHomAltGLsOrGPs_TransversionsOnly
        allCategories$derivedAllelesAll_TV_frac <- allCategories$derivedAllelesAll_TV / allCategories$callableSites_AllCDSSites 
        allCategories$derivedAllelesAll_TV_meanScaled <- allCategories$derivedAllelesAll_TV_frac * allCategories$meanCallableSitesAllSitesAcrossInds # don't need to do factor of 2 bc cancels out
        
        allCategories$HomAlt_meanScaled_frac <- (allCategories$sumHomAltGLsOrGPs / allCategories$callableSites_AllCDSSites)
        allCategories$HomAlt_meanScaled <- (allCategories$sumHomAltGLsOrGPs / allCategories$callableSites_AllCDSSites)*allCategories$meanCallableSitesAllSitesAcrossInds
        
        allCategories$HomAlt_meanScaled_TV_frac <- (allCategories$sumHomAltGLsOrGPs_TransversionsOnly / allCategories$callableSites_AllCDSSites)
        allCategories$HomAlt_meanScaled_TV <- (allCategories$sumHomAltGLsOrGPs_TransversionsOnly / allCategories$callableSites_AllCDSSites)*allCategories$meanCallableSitesAllSitesAcrossInds
        # write out per set of filters:
        write.table(allCategories,paste(data.dir,date,"/angsdOut.mappedTo",ref,".hetHomTotals.",type,".ProbCutoff.",ProbCutoff,".DepthCutoff.",DepthCutoff,".minInd.",minInd,".",date,".AllPointEstimates.allCategories.txt",sep=""),row.names = F,col.names=T,quote=F,sep="\t")
        allInputsAllFilters <- rbind(allInputsAllFilters,allCategories)
        
      }}}
  write.table(allInputsAllFilters,paste(data.dir,date,"/angsdOut.mappedTo",ref,".hetHomTotals.",type,".AllPointEstimates.allCategories.ConcatAllDiffFilterOptions.txt",sep=""),sep="\t",row.names = F,col.names=T,quote=F)
}
# Get overall callable sites:
#allInputs$derivedAllelesAll <- allInputs$sumHetGLsOrGPs + 2*allInputs$sumHomAltGLsOrGPs
# allInputs$derivedAllelesAll_frac <- allInputs$derivedAllelesAll / allInputs$callableSites_AllCDSSites
# allInputs$derivedAllelesAll_meanScaled <- allInputs$derivedAllelesAll_frac * allInputs$meanCallableSitesAllSitesAcrossInds # don't need to do factor of 2 bc cancels out
# # transversions:
# allInputs$derivedAllelesAll_TV <- allInputs$sumHetGLsOrGPs_TransversionsOnly + 2*allInputs$sumHomAltGLsOrGPs_TransversionsOnly
# allInputs$derivedAllelesAll_TV_frac <- allInputs$derivedAllelesAll_TV / allInputs$callableSites_AllCDSSites 
# allInputs$derivedAllelesAll_TV_meanScaled <- allInputs$derivedAllelesAll_TV_frac * allInputs$meanCallableSitesAllSitesAcrossInds # don't need to do factor of 2 bc cancels out
# 
# allInputs$HomAlt_meanScaled_frac <- (allInputs$sumHomAltGLsOrGPs / allInputs$callableSites_AllCDSSites)
# allInputs$HomAlt_meanScaled <- (allInputs$sumHomAltGLsOrGPs / allInputs$callableSites_AllCDSSites)*allInputs$meanCallableSitesAllSitesAcrossInds
# 
# allInputs$HomAlt_meanScaled_TV_frac <- (allInputs$sumHomAltGLsOrGPs_TransversionsOnly / allInputs$callableSites_AllCDSSites)
# allInputs$HomAlt_meanScaled_TV <- (allInputs$sumHomAltGLsOrGPs_TransversionsOnly / allInputs$callableSites_AllCDSSites)*allInputs$meanCallableSitesAllSitesAcrossInds


######### plot point ests ###########
ggplot(allInputs,aes(x=ancOrModern,y=derivedAllelesAll_TV_meanScaled))+
  #geom_violin()+
  geom_point()+
  theme_bw()+
  facet_grid(category~Filter_PerIndividualDepthMinimum~date)+
  ggtitle("Derived Alleles -- scaled by mean sites. Transversions Only!\nvalues are low a) because it's transversions only; b) because avg number of called sites is low due to ancient samples")

### get total homozygous genotypes:

ggplot(allInputsAllFilters,aes(x=ancOrModern,y=HomAlt_meanScaled_TV_frac))+
  #geom_boxplot()+
  geom_point()+
  theme_bw()+
  facet_grid(category~Filter_PerIndividualDepthMinimum~date,scales="free")+
  ggtitle("Transversions only! Fraction of homozygous genotypes")

######### need to do bootstraps of some sort ##################