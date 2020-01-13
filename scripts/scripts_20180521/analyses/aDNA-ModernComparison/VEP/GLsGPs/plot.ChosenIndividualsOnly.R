#### want to plot the aDNA with only one individual ##
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/VEP/GLsGPs/compareMisSynDists_withBootstraps/OnlyTwoInds/"
plot.dir=paste(data.dir,"plots/")
dir.create(plot.dir)
date="20190701"
covs=c("highcov","lowcov")
categories=c("synonymous","missense","stopgained","CDS")
depths=c(1,2,4)
######## these inputs just took the highest cov aDNA sample and one modern CA sample and only looks at sites where they are both called at different filters ####### 
allInputs <- data.frame()
for(cov in covs){
  for(cat in categories){
    for(depth in depths){
      input <- read.table(paste(data.dir,date,"-",cov,"-AFprior-MajorMinor4/angsdOut.mappedTomfur.hetHomTotals.GPs.ProbCutoff.0.95.DepthCutoff.",depth,".minInd.ChosenIndsOnly.",date,"-",cov,"-AFprior-MajorMinor4.",cat,".txt",sep=""),header=T,sep="\t")
      input$coverage <- cov
      input$category <- cat
      input$derivedAllelesTV <- (2*input$sumHomAltGLsOrGPs_TransversionsOnly)+input$sumHetGLsOrGPs_TransversionsOnly
      input$label <- NA
      input[input$sample=="A30_Elut_CA_SM_35_SN1_CAP",]$label <- "Ancient"
      input[input$sample %in% c("116_Elut_CA_307_downsamp","116_Elut_CA_307"),]$label <- "Modern"
      # plot (TV only): derived alelles (TV), hom alt, het between ancient and modern
      allInputs <- rbind(allInputs,input)
    }}}
head(allInputs)

##### order factors:
allInputs$category <- factor(allInputs$category,levels=c("CDS","synonymous","missense","stopgained"))
# abbreviated sample names:
allInputs$shortSamplename <- unlist(lapply(strsplit(as.character(allInputs$sample),"_"),"[",1))
allInputs$shortSamplename <- factor(allInputs$shortSamplename,levels=c("A30","116"))
### plot for different depths and coverages :
for(cov in covs){
  for(depth in depths){
    toPlot <- allInputs[allInputs$coverage==cov & allInputs$Filter_PerIndividualDepthMinimum==depth,]
    # derived alleles #
    plot1 <- ggplot(toPlot,aes(x=shortSamplename,y=derivedAllelesTV,color=label))+
      geom_point()+
      xlab("")+
      facet_wrap(~category,scales="free")+
      theme_bw()+
      ggtitle(paste("Derived Alleles (TV Only)\n",cov,"; MinDepthPerInd=",depth,sep=""))+
      theme(legend.position = "none",text = element_text(size=14))
    ggsave(paste(plot.dir,cov,".minDP",depth,".derivedAlleles.TVOnly.pdf",sep=""),plot1,height=5,width=4)
    # hom-alt #
    plot2 <- ggplot(toPlot,aes(x=shortSamplename,y=sumHomAltGLsOrGPs_TransversionsOnly,color=label))+
      geom_point()+
      xlab("")+
      facet_wrap(~category,scales="free")+
      theme_bw()+
      ggtitle(paste("Hom-Alt Genotypes (TV Only)\n",cov,"; MinDepthPerInd=",depth,sep=""))+
      theme(legend.position = "none")
    ggsave(paste(plot.dir,cov,".minDP",depth,".homAltGTs.TVOnly.pdf",sep=""),plot2,height=5,width=4)
    plot3 <- ggplot(toPlot,aes(x=shortSamplename,y=sumHetGLsOrGPs_TransversionsOnly,color=label))+
      geom_point()+
      xlab("")+
      facet_wrap(~category,scales="free")+
      theme_bw()+
      ggtitle(paste("Het Genotypes (TV Only)\n",cov,"; MinDepthPerInd=",depth,sep=""))+
      theme(legend.position = "none")
    ggsave(paste(plot.dir,cov,".minDP",depth,".hetGTs.TVOnly.pdf",sep=""),plot3,height=5,width=4)

}}


### also want to plot total sites for each category 

plot4 <- ggplot(allInputs[allInputs$category=="CDS",],aes(x=Filter_PerIndividualDepthMinimum,y=callableSites))+
  geom_col()+
  theme_bw()+
  facet_wrap(~coverage)+
  ggtitle("Total CDS sites covered by both individuals")
plot4

ggsave(paste(plot.dir,"jointlycalledsites.pdf"),plot4,height=5,width=7)







