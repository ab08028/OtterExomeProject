####### comparing depth of 1 vs 2 with GP prob of 0.95 only and mfur only NEUTRAL ONLY (bc fast) ##### 
require(ggplot2)
require(reshape2)
#dates=c("20190511","20190513-highcov","20190513-lowcov","20190513-highcov-minInd","20190513-lowcov-minInd")
dates=c("20190524-highcov","20190524-lowcov")
priors="AFprior"
overallDir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/Heterozygosity/"
probCutoffs=0.95
depthCutoffs=c(1,2)
minInds=c(2) # also have 0 somewhere, but for now focus on 2, then look at zero to look at influence of that
refs="mfur"
allInputs=data.frame()

for(angsdDate in dates){
    angsdDatePrior=paste(angsdDate,prior,sep="-")
    data.dir=paste(overallDir,angsdDatePrior,"/",sep="")
    for(ref in refs){
      for(probCutoff in probCutoffs){
        for(depthCutoff in depthCutoffs){
          for(minInd in minInds){
        input <- read.table(paste(data.dir,"angsdOut.mappedTo",ref,".superfile.GPs.neutralOnly.hetHomTotals.ProbCutoff.",probCutoff,".DepthCutoff.",depthCutoff,".minInd.",minInd,".",angsdDatePrior,".txt",sep=""),header=T,sep="\t")
        input$date <- angsdDate
        input$prior <- prior
        input$label <- "modern"
        input[grep("^A",input$sample),]$label <- "ancient"
        # check if any of the sample IDs contain "downsamp":
        if(any(grepl("downsamp",input$sample))){
          input[grep("downsamp",input$sample),]$label <- "modern-downsampled"
        }
        input$reference <- ref
        allInputs = rbind(allInputs,input)
      }
    }
  }
}}

input_melt <- melt(allInputs,measure.vars = c("HetPerSite","HetPerSite_TransversionsOnly"),id.vars = c("sample","label","reference","date","prior","Filter_PerIndividualDepthMinimum","Filter_ProbThresholdForCallableSite","Filter_minIndsPerSite"))

# label transitions+transv vs transv only
input_melt$label2 <- "Transitions+Transversions"
input_melt[input_melt$variable=="HetPerSite_TransversionsOnly",]$label2 <- "TransversionsOnly"



########### okay goal of this plot is to see if the callable sites still scales linearly with the depth =2 ###########
callableSitesVsHet <- allInputs[,c("sample","callableSites","HetPerSite","HetPerSite_TransversionsOnly","Filter_PerIndividualDepthMinimum","Filter_minIndsPerSite","label","reference","date")]
p0a <- ggplot(callableSitesVsHet,aes(x=callableSites,y=HetPerSite,color=label))+
  geom_point()+
  facet_wrap(date~Filter_PerIndividualDepthMinimum)+
  theme_bw()+
  ggtitle("Impact of minimum req'd read depth (1 or 2 reads)\nNeutralOnly")
p0a
p0b <- ggplot(callableSitesVsHet,aes(x=callableSites,y=HetPerSite_TransversionsOnly,color=label))+
  geom_point()+
  facet_wrap(date~Filter_PerIndividualDepthMinimum)+
  theme_bw()+
  ggtitle("Impact of minimum req'd read depth (1 or 2 reads)\nTransversionsOnly\nNeutralOnly")
p0b 

p1 <- ggplot(input_melt[input_melt$variable=="HetPerSite_TransversionsOnly",],aes(x=label,y=value,color=as.factor(Filter_PerIndividualDepthMinimum)))+
  geom_violin(position = position_dodge(0.5))+
  geom_point(position=position_dodge(0.5))+
  facet_wrap(~date,scales="free_x")+
  theme_bw()
p1
