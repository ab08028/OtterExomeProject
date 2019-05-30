### R plot pis from parsing Beagle:
# this is only based on 200K sites
require(ggplot2)
require(reshape2)
#dates=c("20190511","20190513-highcov","20190513-lowcov","20190513-highcov-minInd","20190513-lowcov-minInd")
dates=c("20190524-highcov","20190524-lowcov","20190524-highcov","20190524-lowcov")
priors=c("AFprior","UNIFprior")
overallDir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/Heterozygosity/"
probCutoffs=c(0.5,0.95)
depthCutoff=1

allInputs=data.frame()

for(angsdDate in dates) {
  for(prior in priors){
    angsdDatePrior=paste(angsdDate,prior,sep="-")
    data.dir=paste(overallDir,angsdDatePrior,"/",sep="")
    for(ref in c("mfur","elut")){
      for(probCutoff in probCutoffs){
        input <- read.table(paste(data.dir,"angsdOut.mappedTo",ref,".hetFromPost.ProbCutoff.",probCutoff,".DepthCutoff.",depthCutoff,".",angsdDatePrior,".txt",sep=""),header=T,sep="\t")
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
}

input_melt <- melt(allInputs,measure.vars = c("HetPerSite","HetPerSite_transversionsOnly"),id.vars = c("sample","label","reference","date","prior","PerIndividualDepthMinimum","ProbThresholdForCallableSite"))

# label transitions+transv vs transv only
input_melt$label2 <- "Transitions+Transversions"
input_melt[input_melt$variable=="HetPerSite_transversionsOnly",]$label2 <- "TransversionsOnly"

p0 <- ggplot(input_melt,aes(x=sample,y=value,fill=label2))+
  geom_bar(stat="identity",position="dodge")+
  coord_flip()+
  theme_bw()+
  facet_grid(reference~ProbThresholdForCallableSite~date~prior)+
  theme(legend.title=element_blank())
p0

ggsave(paste(overallDir,"perIndividualHets.allInds.pdf"),device="pdf",height=5,width=9)
# Based on this, UNIF distribution way overestimates things. Don't use. Stick to AF, but maybe use a more perind filter.
####### AF prior only: 
p1 <- ggplot(input_melt[input_melt$prior=="AFprior",],aes(x=label,y=value,fill=label))+
  geom_violin(position=position_dodge(.5))+
  geom_point(position=position_dodge(.5),size = 1)+
  theme_bw()+
  ggtitle("Comparing Heterozygosity")+
  theme(legend.title = element_blank(),axis.text = element_text(size=14),legend.text = element_text(size=14),legend.background = element_rect(fill = "transparent"))+
  xlab("")+
  facet_grid(ProbThresholdForCallableSite~reference~label2~date,scales="free")
p1
ggsave(paste(overallDir,"ancientVsModernHets.allInds.pdf"),p1,device="pdf",height=5,width=11)

# this plot sucks: 
p2 <- ggplot(input_melt[input_melt$variable=="HetPerSite_transversionsOnly",],aes(x=label,y=value,fill=interaction(label,date)))+
  geom_violin(position=position_dodge(.5))+
  geom_point(position=position_dodge(.5),size = 1)+
  theme_bw()+
  ggtitle("Comparing Heterozygosity -- Transversions Only")+
  theme(legend.title = element_blank(),axis.text = element_text(size=14))+
  xlab("")+
  facet_wrap(~reference)
p2
ggsave(paste(overallDir,"ancientVsModernHets.allInds.TransvesionsOnly.CompareCoverage.pdf"),p2,device="pdf",height=5,width=11)


########## you are here ## plots don't look great -- try to improve them !! ######
# things to separate by:
# transitions + transversions / transversions
# elut vs mfur
# condition (what priors were based on)

# label whether dataset contains high cov or low cov or both:
allInputs$label2 <- NA
allInputs[allInputs$date=="20190511",]$label2 <- "AllInds"
allInputs[allInputs$date=="20190513-lowcov-minInd",]$label2 <- "LowCov"
allInputs[allInputs$date=="20190513-lowcov",]$label2 <- "LowCov"
allInputs[allInputs$date=="20190513-highcov",]$label2 <- "HighCov"
allInputs[allInputs$date=="20190513-highcov-minInd",]$label2 <- "HighCov"

# label whether dataset had minInd filter
allInputs$label3 <- NA
allInputs[allInputs$date=="20190511",]$label3 <- "NoMindInd"
allInputs[allInputs$date=="20190513-lowcov-minInd",]$label3 <- "MinInd5"
allInputs[allInputs$date=="20190513-lowcov",]$label3 <- "NoMindInd"
allInputs[allInputs$date=="20190513-highcov",]$label3 <- "NoMindInd"
allInputs[allInputs$date=="20190513-highcov-minInd",]$label3 <- "MinInd5"

# tranisitons+transversions:
hetPlot1a <- ggplot(allInputs,aes(x=label,y=HetPerSite,color=label3,shape=label2))+
  geom_point(size=4,alpha=0.5,position=position_dodge(width=0.4))+
  #geom_violin(aes(x=date,y=HetPerSite_transversionsOnly,fill=label),alpha=0.5)+
  facet_wrap(~reference)+
  theme_bw()+
  ggtitle("Transtions+Transversions")
hetPlot1a
ggsave(paste(overallDir,"InfluenceOfSettingsOnHet.pdf",sep=""),hetPlot1a,width=7,height=7)

# transversions only:
hetPlot1b <- ggplot(allInputs,aes(x=label,y=HetPerSite_transversionsOnly,color=label3,shape=label2))+
  geom_point(size=4,alpha=0.5,position=position_dodge(width=0.4))+
  #geom_violin(aes(x=date,y=HetPerSite_transversionsOnly,fill=label),alpha=0.5)+
  facet_wrap(~reference)+
  theme_bw()+
  ggtitle("Transversions Only")
hetPlot1b
ggsave(paste(overallDir,"InfluenceOfSettingsOnHet.TransversionsOnly.pdf",sep=""),hetPlot1b,width=7,height=7)


require(scales)
# plot callable sites
callableSitesPlot <- ggplot(allInputs,aes(x=label2,y=callableSites,fill=label3))+
  geom_bar(stat="identity",position="dodge")+
  #coord_flip()+
  theme_bw()+
  theme(legend.title=element_blank())+
  scale_y_continuous()+
  xlab("Dataset")+
  ggtitle("Comparison of number of passing sites depending on individuals in dataset")
callableSitesPlot
ggsave(paste(overallDir,"CallableSitesDependingOnDataset.pdf",sep=""),callableSitesPlot,width=7,height=7)
