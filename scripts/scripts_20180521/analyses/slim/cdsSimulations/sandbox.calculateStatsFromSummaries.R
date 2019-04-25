####### sandbox #######
# play with summaries 
require(dplyr)
require(ggplot2)
numChunk=20 # chunks per replicate
numReps=20 # number of reps
sim.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/concattedSummaries/"
plot.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/slim/"
# ignore for now: "COM/1D.3Epoch.1.5Mb.cds/20190404/h_0/","AK/1D.2Epoch.1.5Mb.cds/20190404/h_0/","AL/1D.2Epoch.1.5Mb.cds/20190404/h_0/"
rundate=20190424 # set of simulations you're interested in (if is across different rundates you can list popsModelsRundates explicitly)
hs=c("0.5","0") # set of hs you're interested in
popMods=c("AL/1D.2Epoch.1.5Mb.cds", "AK/1D.2Epoch.1.5Mb.cds") # population and corresponding models you're interested in
popmodels=list()
for(i in popMods){
  for(h in hs){
    pm=paste(i,"/",rundate,"/h_",h,"/",sep="")
    popmodels=c(popmodels,pm)
  }
}

#popmodels=c("AK/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/", "AL/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/", "CA/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/", "KUR/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/")
# weird bottleneck model -- not sure what to make of it "COM/1D.3Epoch.1.5Mb.cds/20190423/h_0.5/"
#numreps=20 # number of reps 
# go through pre and post 
allSummaries = data.frame()
# this will work once bot pre and post have pop size column
for(pm in popmodels){ # go through each population-model dir
  for(state in c("Post","Pre")){  # do separtely for pre and post contraction 
    #for(rep in seq(1,numreps)){ # don't go through each replicate because rep numbers differ between populations/models now that I select the first passing 20 replicates
    # NOTE: each pop/model will have e.g. 20 replicates, but they may be diff numbers
    # This is because I do 25 replicates, and pick out the first 20 to finish (to account for random Hoffman failures). Each replicate will be fully complete with the requisite number of chunks, so don't have to worry about total number of sites simulated etc. that will be fully complete within passing replicates
    replicates = list.files(paste(sim.dir,pm,sep = ""),pattern=paste("replicate_.*",state,sep=""),full.names = T)
    for(infile in replicates) {
      input <- read.table(infile,header=T,sep=",")
      # make sure there are 20 chunks and 20 reps:
      checkChunk=length(input$chunk) # should be numChunk (eg 20) chunks per replicate
      if(checkChunk!=numChunk){
        print(paste("NOT ENOUGH CHUNKS IN THIS REPLICATE: ",infile))
        break
      }
      input$state <- state
      input$popModel <- pm
      input$population <- unlist(lapply(strsplit(input$popModel,"/"),"[",1))
      input$model <- unlist(lapply(strsplit(input$popModel,"/"),"[",2))
      input$modelDate <- unlist(lapply(strsplit(input$popModel,"/"),"[",3))
      input$dominance_h <- unlist(lapply(strsplit(input$popModel,"/"),"[",4))
      input$state <- paste("simulated: ",state,"-Contraction",sep="")
      input$category <- "simulated"
      allSummaries = rbind(input,allSummaries)
    }
  }
}


allSummaries$mutLabel <- NA
allSummaries[allSummaries$type=="m1",]$mutLabel <- "synonymous"
allSummaries[allSummaries$type=="m2",]$mutLabel <- "missense"

############# add s categories ##########
allSummaries$sCat <- NA
# JAR categories from her 2018 paper
allSummaries[allSummaries$s >= -1 & allSummaries$s < -0.01,]$sCat <- "strongly deleterious"
allSummaries[allSummaries$s >= -0.01 & allSummaries$s < -0.001,]$sCat <- "moderately deleterious"
allSummaries[allSummaries$s >= -0.001 & allSummaries$s < 0,]$sCat <- "weakly deleterious"
allSummaries[allSummaries$s==0,]$sCat <- "neutral"
#strongly deleterious, −1 ≤ s <−0.01; ; (C) moderately deleterious, −0.01 ≤ s <−0.001; (D) weakly deleterious, −0.001 ≤ s < 0; and (E) neutral, s = 0. S

###### want to get the mean number of alleles per category per individual
# it's a population wide-mean so you divide by the popsize at the time you sampled the population (which is a column in the output)
# pop model includes dominance_h

# group by replicate?
totalDerivedAlleles <- allSummaries %>%
  group_by(replicate,popModel,population,state,sCat,model,dominance_h,popsizeDIP) %>%
  tally((sum(p1numhom)+2*sum(p1numhet))) # this works to get total 


######## DO NOT DO THE DIVIDING BY POP SIZE IN THIS DPLYR STEP ^; for some reason really messes up estimates (unsure why right now); do it in the following step:  ########
totalDerivedAlleles$meanDerivedAlleles <- totalDerivedAlleles$n / totalDerivedAlleles$popsizeDIP


###### then get the mean across replicates ########
totalDerivedAlleles$state <- factor(totalDerivedAlleles$state,levels=c("simulated: Pre-Contraction","simulated: Post-Contraction"))
totalDerivedAlleles$dominance_h <- factor(totalDerivedAlleles$dominance_h,levels=c("h_0.5","h_0"))
totalDerivedAlleles$sCat <- factor(totalDerivedAlleles$sCat,levels=c("neutral","weakly deleterious","moderately deleterious","strongly deleterious"))


#p1 <- ggplot(totalDerivedAlleles_withMeans,aes(x=population,y=MeanAcrossReps,fill=state))+
#  geom_bar(stat="identity",position="dodge")+
#  facet_wrap(~sCat~dominance_h~model,scales="free")+
#  theme_bw()

#p1
#ggsave(paste(plot.dir,"meanDerivedAllelesPerIndividual.Jacqueline-like.pdf",sep=""),p1,device="pdf",height=7,width=9)


### try boxplot

# ggplot(totalDerivedAlleles,aes(x=population,y=meanDerivedAlleles,fill=state))+
#   #geom_point(stat="identity",position=position_dodge(width=0.9),shape=5) +
#   geom_boxplot()+
#   #geom_point(position=position_dodge(width=0.9),data=totalDerivedAlleles,aes(x=population,y=meanDerivedAlleles,color=state))+
#   facet_wrap(~sCat~dominance_h~model,scales="free")


## This is wrong! this isnt' mean! is stacking geoms 
repMeans <- totalDerivedAlleles %>%
  group_by(population,model,popModel,sCat,dominance_h,popsizeDIP,state) %>%
  summarise(meanAcrossReps=mean(meanDerivedAlleles))

p1 <- ggplot(totalDerivedAlleles,aes(x=population,y=meanDerivedAlleles))+
  # points for individual runs
  geom_point(position=position_dodge(width=0.9),aes(color=state),shape=1,alpha=0.8)+
  # bar for mean
  geom_bar(data=repMeans,aes(x=population,y=meanAcrossReps,fill=state),stat="identity",position=position_dodge(width=0.9),alpha=0.8) +
  scale_color_manual(values=c("black","black"))+  
  scale_fill_manual(values=c("darkblue","darkred"))+
  facet_wrap(~model~sCat~dominance_h,scales="free")+
  theme_bw()+
  ggtitle(paste("bar is mean of population means across ", numReps," replicates\nDots are individual replicate population means"))
p1
ggsave(paste(plot.dir,"meanDerivedAllelesPerIndividual.acrossReps.pdf",sep=""),p1,device="pdf",height=8,width=10)
