######### Look at polymorphic sites per population #############
require(ggplot2)

data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/countingPolymorphicSites/"
dates=c("20190612")
datasets=c("highcov") # low cov is still running ,"lowcov")
MaxMissing=c(0,1) # 0 is more conservative; 1 is less
refs=c("mfur","elut")
all.inputs=data.frame()
for(date in dates){
  for(dataset in datasets){
  for(ref in refs){
    for(missing in MaxMissing){
      inputfile <- paste(data.dir,date,"-",dataset,"-pseudoHaps/",dataset,".PolymorphicCounts.mappedTo",ref,".maxMiss.",missing,".txt",sep="")
      input <- read.table(inputfile,header=T,sep="\t")
      input$ref <- ref
      input$dataset <- dataset
      input$date <- date
      input$MaxMissingInd <- as.factor(input$MaxMissingInd)
      all.inputs <- rbind(all.inputs,input)
    }}}}
all.inputs$polyTVFraction <- all.inputs$PolymorphicSitesTV / all.inputs$sitesCalled
all.inputs$polyAllFraction <- all.inputs$PolymorphicSitesAll / all.inputs$sitesCalled

# oder the factors:
all.inputs$pop <- factor(all.inputs$pop,levels=c("ancient","AK","CA"))
### get the fraction of sites that are polymorphic 
p1 <- ggplot(all.inputs,aes(x=pop,y=polyTVFraction,fill=MaxMissingInd))+
  geom_bar(stat="identity",position="dodge",alpha=0.7)+
  facet_wrap(~ref~dataset)+
  theme_bw()+
  ggtitle("Fraction of called sites that are polymorphic: Transversions Only\nPseudohaploids")
p1
ggsave(paste(data.dir,"fracPolymorphicsites.TransversionOnly.pseudohaps.anc.ca.ak.pdf",sep=""),p1,height=5,width=7)

# total called:
ggplot(all.inputs, aes(x=pop,y=sitesCalled,fill=MaxMissingInd))+
  geom_bar(stat="identity",position="dodge",alpha=0.7)+
  facet_wrap(~ref~dataset)+
  theme_bw()
  