####### sandbox #######
# play with summaries 
require(dplyr)
require(ggplot2)
sim.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/concattedSummaries/"
plot.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/slim"
# ignore for now: "COM/1D.3Epoch.1.5Mb.cds/20190404/h_0/","AK/1D.2Epoch.1.5Mb.cds/20190404/h_0/","AL/1D.2Epoch.1.5Mb.cds/20190404/h_0/"
popmodels=c("AK/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/", "AL/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/", "CA/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/", "KUR/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/")
# weird bottleneck model -- not sure what to make of it "COM/1D.3Epoch.1.5Mb.cds/20190423/h_0.5/"
numreps=1 # number of reps 
# go through pre and post 
allSummaries = data.frame()
# this will work once bot pre and post have pop size column
for(pm in popmodels){ # go through each population-model dir
  for(state in c("Post","Pre")){  # do separtely for pre and post contraction 
    for(rep in seq(1,numreps)){ # go through each replicate
      input <- read.table(paste(sim.dir,pm,"rep.",rep,".slim.output.",state,"Contraction.allConcatted.summary.txt.gz",sep=""),sep=",",header=T)
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


totalDerivedAlleles <- allSummaries %>%
  group_by(popModel,population,state,sCat,model,dominance_h,popsizeDIP) %>%
  tally((sum(p1numhom)+2*sum(p1numhet))) # this works to get total 
######## DO NOT DO THE DIVIDING BY POP SIZE IN THIS DPLYR STEP ^; for some reason really messes up estimates (unsure why right now); do it in the following step:  ########
totalDerivedAlleles$meanDerivedAlleles <- totalDerivedAlleles$n / totalDerivedAlleles$popsizeDIP

# order factors for plotting:
totalDerivedAlleles$state <- factor(totalDerivedAlleles$state,levels=c("simulated: Pre-Contraction","simulated: Post-Contraction"))
totalDerivedAlleles$sCat <- factor(totalDerivedAlleles$sCat,levels=c("neutral","weakly deleterious","moderately deleterious","strongly deleterious"))
p1 <- ggplot(totalDerivedAlleles,aes(x=interaction(population),y=meanDerivedAlleles,fill=state))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~sCat~dominance_h~model,scales="free")+
  theme_bw()+
  ggtitle("note: this is a single simulation replicate\nand the CA and KUR models are dubious (crash to ~6 for 1 gen which isn't v. plausible)")

p1

