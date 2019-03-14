require(ggplot2)
require(gridExtra)
####### sandbox trying to plot dadi results 
slimmodel="1D.2Epoch" # model you simulated under
slimdate=20190125 # date you ran slim
pop="generic" # what you named population in slim simulation
Nanc_true=4000 # nanc you simulated
nu_true=30 # nu you simulated
T_true=10 # gen you simulated
# need to show a couple things
# 1) whether 1Epoch or 2Epoch fits better for the same replicate
# 2) if parameters converge for 2Epoch
# 3) what the range of those parameter estimates are 
dadimodel="1D.2Epoch" # model you inferred in dadi
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/neutralSimulations/",slimmodel,"/",slimdate,"/dadiInfBasedOnSlim/allDadiResultsConcatted/dadiInfModel_",dadimodel,"/",sep = "")
data.dir

# gather up best LLs from each replicate
allbest <- NULL
for(rep in seq(1,11)){
  results <- read.table(paste(data.dir,pop,".rep.",rep,".dadi.inf.",dadimodel,".all.output.concatted.txt",sep=""),sep = "\t",header=T,stringsAsFactors = F)
  # select best LL
  best <- results[results$LL==max(results$LL),]
  best$rep=as.character(rep)
  # add to allbest:
  allbest <- rbind(allbest,best)
}

allresults <- NULL
for(rep in seq(1,11)){
  results <- read.table(paste(data.dir,pop,".rep.",rep,".dadi.inf.",dadimodel,".all.output.concatted.txt",sep=""),sep = "\t",header=T,stringsAsFactors = F)
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

ggplot(allbest,aes(x=nu_scaled_dip))+
  geom_histogram(alpha=0.6)+
  geom_vline(xintercept = nu_true,color="darkred")+
  ggtitle("Best (out of 50) LL estimate of nu\nAcross simulation replicates\n(dadi 1D.2Epoch inferred from Slim 1D.2Epoch")+
  theme_bw()

ggplot(allbest,aes(x=T_scaled_gen))+
  geom_histogram(alpha=0.6)+
  geom_vline(xintercept = T_true,color="darkred")+
  ggtitle("Best (out of 50) LL estimate of nu\nAcross simulation replicates\n(dadi 1D.2Epoch inferred from Slim 1D.2Epoch")+
  theme_bw()

################  plot all 2 Epoch results #############

p1 <- ggplot(allresults,aes(x=reorder(rep, sort(as.numeric(rep))),y=Nanc_FromTheta_scaled_dip,color=LL))+
  geom_point()+
  geom_hline(yintercept = Nanc_true)+
  scale_color_gradientn(colours=rainbow(5))+
  ylab("Nanc\n(ancestral size)")+
  xlab("Replicate")
p1
p2 <- ggplot(allresults,aes(x=reorder(rep, sort(as.numeric(rep))),y=T_scaled_gen,color=LL))+
  geom_point()+
  geom_hline(yintercept = T_true)+
  scale_color_gradientn(colours=rainbow(6))+
  scale_y_log10()+
  ylab("T\n(time of contraction)")+
  xlab("Replicate")
p2 
p3 <- ggplot(allresults,aes(x=reorder(rep, sort(as.numeric(rep))),y=nu_scaled_dip,color=LL))+
  geom_point()+
  geom_hline(yintercept = nu_true)+
  scale_color_gradientn(colours=rainbow(6))+
  scale_y_log10()+
  ylab("nu\n(contraction size)")+
  xlab("Replicate")
p3

allP <- grid.arrange(p1,p3,p2)


ggsave(paste(data.dir,"comparingDadiReplicates.pdf",sep=""),allP,width=9,height=6)


################## without log scaling ##########
p2b <- ggplot(allresults,aes(x=reorder(rep, sort(as.numeric(rep))),y=T_scaled_gen,color=LL))+
  geom_point()+
  geom_hline(yintercept = T_true)+
  scale_color_gradientn(colours=rainbow(6))+
  #scale_y_log10()+
  ylab("T\n(time of contraction)")+
  xlab("Replicate")
p2b 
p3b <- ggplot(allresults,aes(x=reorder(rep, sort(as.numeric(rep))),y=nu_scaled_dip,color=LL))+
  geom_point()+
  geom_hline(yintercept = nu_true)+
  scale_color_gradientn(colours=rainbow(6))+
  #scale_y_log10()+
  ylab("nu\n(contraction size)")+
  xlab("Replicate")
p3b

allPb <- grid.arrange(p1,p3b,p2b)


ggsave(paste(data.dir,"comparingDadiReplicates.noLogScale.pdf",sep=""),allPb,width=9,height=6)

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
  geom_point(size=2)+
  geom_point(data=oneResultPer_rep_1Epoch,aes(x=rep,y=LL,color="1Epoch"),size=2)+
  ggtitle("Compare LLs between 1Epoch and 2Epoch Dadi Models")+
  theme(legend.title = element_blank())+xlab("Replicate")

p4
##### gather them all together
#grid.arrange(p1,p3,p2,p4,ncol=1)
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/neutralSimulations/",slimmodel,"/",slimdate,"/dadiInfBasedOnSlim/allDadiResultsConcatted/","comparing.1Epoch.2Epoch.DadiLLs.perReplicate.pdf",sep=""),p4,width=8.5,height=3)
