require(ggplot2)
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/sumGPsGLsPerVEPCategory/"

dates=c("20190701-highcov-AFprior-MajorMinor4","20190701-lowcov-AFprior-MajorMinor4")
minInds=c(1)
ProbCutoffs=c(0.95)
DepthCutoffs=c(1,2,4)
# cds is all coding; missense; synonymous (skipping SG for now since it's a bit weird after using --pick.)
categories=c("missense","synonymous") # something's a bit odd with stopgained, too few? or too many in otter paper perhaps. Maybe --pick reclassifies them, that is most likely. then that's okay. because then it's a parameter choice and it is the same between each otter in genome paper, so even if absolute # is disputable, the *pattern* is the same. Okay. cool.

ref="mfur"
type="GPs"
allInputs <- data.frame()
for(date in dates){
    for(minInd in minInds){
    for(ProbCutoff in ProbCutoffs){
      for(DepthCutoff in DepthCutoffs){
        # get total cds sites:
        CDSinfile=paste(data.dir,date,"/angsdOut.mappedTo",ref,".hetHomTotals.",type,".ProbCutoff.",ProbCutoff,".DepthCutoff.",DepthCutoff,".minInd.",minInd,".",date,".CDS.txt",sep="")
        cdsSites <- read.table(CDSinfile,header=T,sep="\t")
        callableSites <- cdsSites[,c("sample","callableSites")]
        meanCallableSitesAcrossInds <- mean(callableSites$callableSites)
        # merge this with each input
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
        allInputs <- rbind(allInputs,comboInput)
      }}}}}
      
# Get overall callable sites:
allInputs$derivedAllelesAll <- allInputs$sumHetGLsOrGPs + 2*allInputs$sumHomAltGLsOrGPs
allInputs$derivedAllelesAll_frac <- allInputs$derivedAllelesAll / allInputs$callableSites_AllCDSSites
allInputs$derivedAllelesAll_meanScaled <- allInputs$derivedAllelesAll_frac * allInputs$meanCallableSitesAllSitesAcrossInds # don't need to do factor of 2 bc cancels out
# transversions:
allInputs$derivedAllelesAll_TV <- allInputs$sumHetGLsOrGPs_TransversionsOnly + 2*allInputs$sumHomAltGLsOrGPs_TransversionsOnly
allInputs$derivedAllelesAll_TV_frac <- allInputs$derivedAllelesAll_TV / allInputs$callableSites_AllCDSSites 
allInputs$derivedAllelesAll_TV_meanScaled <- allInputs$derivedAllelesAll_TV_frac * allInputs$meanCallableSitesAllSitesAcrossInds # don't need to do factor of 2 bc cancels out

allInputs$HomAlt_meanScaled_frac <- (allInputs$sumHomAltGLsOrGPs / allInputs$callableSites_AllCDSSites)
allInputs$HomAlt_meanScaled <- (allInputs$sumHomAltGLsOrGPs / allInputs$callableSites_AllCDSSites)*allInputs$meanCallableSitesAllSitesAcrossInds

allInputs$HomAlt_meanScaled_TV_frac <- (allInputs$sumHomAltGLsOrGPs_TransversionsOnly / allInputs$callableSites_AllCDSSites)
allInputs$HomAlt_meanScaled_TV <- (allInputs$sumHomAltGLsOrGPs_TransversionsOnly / allInputs$callableSites_AllCDSSites)*allInputs$meanCallableSitesAllSitesAcrossInds

ggplot(allInputs,aes(x=ancOrModern,y=derivedAllelesAll_TV_meanScaled))+
  #geom_violin()+
  geom_point()+
  theme_bw()+
  facet_grid(category~Filter_PerIndividualDepthMinimum~date)+
  ggtitle("Derived Alleles -- scaled by mean sites. Transversions Only!\nvalues are low a) because it's transversions only; b) because avg number of called sites is low due to ancient samples")

### get total homozygous genotypes:

ggplot(allInputs,aes(x=ancOrModern,y=HomAlt_meanScaled_TV_frac))+
  #geom_boxplot()+
  geom_point()+
  theme_bw()+
  facet_grid(category~Filter_PerIndividualDepthMinimum~date,scales="free")+
  ggtitle("Transversions only! Fraction of homozygous genotypes")

######### need to do bootstraps of some sort ##################