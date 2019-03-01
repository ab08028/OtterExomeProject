require(ggplot2)
require(gridExtra)
####### sandbox trying to plot dadi results 
slimmodel="1D.1Epoch" # model you simulated under
slimdate=20190213 # date you ran slim
pop="generic" # what you named population in slim simulation
Nanc_true=4000 # nanc you simulated
# need to show a couple things
# 1) whether 1Epoch or 2Epoch fits better for the same replicate
# 2) if parameters converge for 2Epoch
# 3) what the range of those parameter estimates are 
dadimodel="1D.1Epoch" # model you inferred in dadi
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/neutralSimulations/",slimmodel,"/",slimdate,"/dadiInfBasedOnSlim/allDadiResultsConcatted/dadiInfModel_",dadimodel,"/",sep = "")
data.dir

# gather up best LLs from each replicate
allbest <- NULL
for(rep in seq(1,11)){
  results <- read.table(paste(data.dir,pop,".rep.",rep,".dadi.inf.",dadimodel,".all.output.concatted.txt",sep=""),sep = "\t",header=T,stringsAsFactors = F,fill=T) # need fill=T because 1 epoch models are missing last 2 columns bc no initial params
  # select best LL
  best <- results[results$LL==max(results$LL),]
  best$rep=as.character(rep)
  # add to allbest:
  allbest <- rbind(allbest,best)
}

allresults <- NULL
for(rep in seq(1,11)){
  results <- read.table(paste(data.dir,pop,".rep.",rep,".dadi.inf.",dadimodel,".all.output.concatted.txt",sep=""),sep = "\t",header=T,stringsAsFactors = F,fill = T) # need fill=T because 1 epoch models are missing last 2 columns bc no initial params
  results$rep=as.character(rep)
  # add to allresults:
  allresults <- rbind(allresults,results)
}
########## try plotting
require(ggplot2)
ggplot(allbest,aes(x=Nanc_FromTheta_scaled_dip))+
  geom_histogram(alpha=0.6)+
  geom_vline(xintercept = Nanc_true,color="darkred")+
  ggtitle("Best (out of 50) LL estimate of Nanc\nAcross simulation replicates\n(dadi 1D.2Epoch inferred from Slim 1D.2Epoch")+
  theme_bw()



################  plot all 2 Epoch results #############

p1 <- ggplot(allresults,aes(x=reorder(rep, sort(as.numeric(rep))),y=Nanc_FromTheta_scaled_dip,color=LL))+
  geom_point()+
  geom_hline(yintercept = Nanc_true)+
  scale_color_gradientn(colours=rainbow(5))+
  ylab("Nanc\n(ancestral size)")+
  xlab("Replicate")
p1


ggsave(paste(data.dir,"comparingDadiReplicates.pdf",sep=""),p1,width=9,height=2)


################## For each pair, plot 1Epoch/2Epoch LLs ###################
dadimodel="1D.1Epoch" # model you inferred in dadi
data.dir2=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/neutralSimulations/",slimmodel,"/",slimdate,"/dadiInfBasedOnSlim/allDadiResultsConcatted/dadiInfModel_",dadimodel,"/",sep = "")
oneResultPer_rep_1Epoch <- NULL
# note: it doesn't matter about 50 replicates of the 1Epoch model; it is determinative
for(rep in seq(1,11)){
  results <- read.table(paste(data.dir2,pop,".rep.",rep,".dadi.inf.",dadimodel,".all.output.concatted.txt",sep=""),sep = "\t",header=T,fill = T) # have to fill in because initial param fields are empty
  results$rep=as.character(rep)
  # add to allresults: just choose the first entry ; didn't need 50 replicates because is determinative
  oneResultPer_rep_1Epoch <- rbind(oneResultPer_rep_1Epoch,results[1,])
}

# pair 1Epoch and 2Epoch by replicate
p4 <- ggplot(allresults,aes(x=reorder(rep, sort(as.numeric(rep))),y=LL,color="2Epoch"))+
  geom_point(size=2,shape=1)+
  geom_point(data=oneResultPer_rep_1Epoch,aes(x=rep,y=LL,color="1Epoch"),size=2)+
  ggtitle("Compare LLs between 1Epoch and 2Epoch Dadi Models")+
  theme(legend.title = element_blank())+xlab("Replicate")

p4
##### gather them all together
#grid.arrange(p1,p3,p2,p4,ncol=1)
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/neutralSimulations/",slimmodel,"/",slimdate,"/dadiInfBasedOnSlim/allDadiResultsConcatted/","comparing.1Epoch.2Epoch.DadiLLs.perReplicate.pdf",sep=""),p4,width=8.5,height=2)
