################### Plotting Pooneh's simulations ################
# want to plot generations after bottleneck
# want to facet by additive/recessive and by population 
require(ggplot2)
require(dplyr)
require(RColorBrewer)
require(viridis) # nice colors
colorPal=RColorBrewer::brewer.pal(n=6,name = "Dark2")
colors=list(CA=colorPal[1],BAJ=colorPal[1],AK=colorPal[2],AL=colorPal[3],COM=colorPal[4],KUR=colorPal[5])
migrationColors=list(NoGeneFlow="grey22","1perGen"="cyan4","5perGen"="#8C96C6","10perGen"="#88419D","25for2Gen"="#6E016B","25perGen"="darkgoldenrod4")
textsize=14

########## I. Pre-Post Contraction + exclude Recovery #################
AK_CA_PrePost <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/poonehSimulations/concattedSummaries/LoadPerGeneration.RemovedPreFurTradeFixedVars.ThroughTime.AllReps.CA.AK.txt",sep="\t",header=T)
# want to exclude things that are post-fur trade (recovery) for one set of plots:
AK_CA_PrePost_noRec <- AK_CA_PrePost[(AK_CA_PrePost$population=="AK" & AK_CA_PrePost$generation<36)|(AK_CA_PrePost$population=="CA" & AK_CA_PrePost$generation<26),]
# this now has no recovery periods
# rest should work as before (check it)
AK_CA_PrePost_noRec$label <- ""
AK_CA_PrePost_noRec[AK_CA_PrePost_noRec$h==0,]$label <- "recessive"
AK_CA_PrePost_noRec[AK_CA_PrePost_noRec$h==0.5,]$label <- "additive"
# gather means:
AK_CA_PrePost_noRec_meansOverSimReplicates = AK_CA_PrePost_noRec %>%
  group_by(model,population,label,generation)  %>%
  summarise(meanLoad=mean(L))

######################### p1: Plot the bottleneck load only ( no recovery ) ################
for(dominance in c("recessive","additive")){
  p1 <- ggplot(AK_CA_PrePost_noRec_meansOverSimReplicates[AK_CA_PrePost_noRec_meansOverSimReplicates$label==dominance,],aes(x=generation,y=meanLoad,color=population))+
    geom_line(data=AK_CA_PrePost_noRec[AK_CA_PrePost_noRec$label==dominance,],aes(group=replicate,x=generation,y=L),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    facet_wrap(~population,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),legend.position = "none",strip.background = element_rect("transparent"))+
    xlab("generations since bottleneck")+
    ylab("mean genetic load")+
    scale_color_manual(values=unlist(colors))+
    ggtitle(paste(dominance," mutations",sep=""))
  p1
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/BOTTLENECK.",dominance,".geneticLoad.AK.CA.RemovedPreFurTradeFixedVars.pdf",sep=""),p1,height=2,width=3)
  
  ############ plot Alaska only ###########
  p1b <- ggplot(AK_CA_PrePost_noRec_meansOverSimReplicates[AK_CA_PrePost_noRec_meansOverSimReplicates$label==dominance & AK_CA_PrePost_noRec_meansOverSimReplicates$population=="AK",],aes(x=generation,y=meanLoad,color=population))+
    geom_line(data=AK_CA_PrePost_noRec[AK_CA_PrePost_noRec$label==dominance & AK_CA_PrePost_noRec$population=="AK",],aes(group=replicate,x=generation,y=L),alpha=0.3,size=0.5)+
    geom_line(size=1.5)+
    theme_bw()+
    facet_wrap(~population,ncol=2,scales="free_x")+
    theme(text = element_text(size=12),legend.position = "none",strip.background = element_rect("transparent"))+
    xlab("generations since bottleneck")+
    ylab("mean genetic load")+
    scale_color_manual(values=unlist(colors))+
    ggtitle(paste(dominance," mutations",sep=""))
  p1b
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/BOTTLENECK.",dominance,".geneticLoad.AK.Only.RemovedPreFurTradeFixedVars.pdf",sep=""),p1b,height=2,width=3)

######## plot with dominance faceted AK only  #########
AK_CA_PrePost_noRec_meansOverSimReplicates$label2 <- paste(AK_CA_PrePost_noRec_meansOverSimReplicates$label," mutations",sep="")
  AK_CA_PrePost_noRec_meansOverSimReplicates$label2 <- factor(AK_CA_PrePost_noRec_meansOverSimReplicates$label2, levels=c("recessive mutations","additive mutations"))
AK_CA_PrePost_noRec$label2 <- paste(AK_CA_PrePost_noRec$label," mutations",sep="")
AK_CA_PrePost_noRec$label2 <- factor(AK_CA_PrePost_noRec$label2, levels=c("recessive mutations","additive mutations"))

p1c <- ggplot(AK_CA_PrePost_noRec_meansOverSimReplicates[AK_CA_PrePost_noRec_meansOverSimReplicates$population=="AK",],aes(x=generation,y=meanLoad,color=population))+
  geom_line(data=AK_CA_PrePost_noRec[AK_CA_PrePost_noRec$population=="AK",],aes(group=replicate,x=generation,y=L),alpha=0.3,size=0.5)+
  geom_line(size=1.5)+
  theme_bw()+
  facet_wrap(~label2,scales="free_y",ncol=1)+
  theme(text = element_text(size=12),legend.position = "none",strip.background = element_rect("transparent"))+
  xlab("generations since bottleneck")+
  ylab("mean genetic load")+
  scale_color_manual(values=unlist(colors))+
  theme(text = element_text(size=11),strip.background = element_blank(),legend.title = element_blank(),strip.text = element_text(hjust=0,size=textsize))
p1c
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/BOTTLENECK.FacetedDominance.geneticLoad.AK.Only.RemovedPreFurTradeFixedVars.pdf",sep=""),p1c,height=4,width=3)
}
########## II. Pre-Post Contraction + Recovery #################
AK_CA_PrePost$label <- ""
AK_CA_PrePost[AK_CA_PrePost$h==0,]$label <- "recessive"
AK_CA_PrePost[AK_CA_PrePost$h==0.5,]$label <- "additive"
# gather means:
AK_CA_Rec_Combo_meansOverSimReplicates = AK_CA_PrePost %>%
  group_by(model,population,label,generation)  %>%
  summarise(meanLoad=mean(L))

######################### p3: Plot with recovery load ################
hLines=data.frame(population=c("AK","AK","CA","CA"),intercept=c(0,36,0,26))
for(dominance in c("recessive","additive")){
  
  p3 <- ggplot(AK_CA_Rec_Combo_meansOverSimReplicates[AK_CA_Rec_Combo_meansOverSimReplicates$label==dominance,],aes(x=generation,y=meanLoad,color=population))+
    geom_line(data=AK_CA_PrePost[AK_CA_PrePost$label==dominance,],aes(group=replicate,x=generation,y=L),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    geom_vline(data=hLines,aes(xintercept = intercept),linetype="dashed",color="black")+
    facet_wrap(~population,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),legend.position = "none",strip.background = element_rect("transparent"))+
    xlab("generations since bottleneck")+
    ylab("mean genetic load")+
    scale_color_manual(values=unlist(colors))+
    ggtitle(paste(dominance," mutations",sep=""))
    
  p3
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/RECOVERY.",dominance,".geneticLoad.AK.CA.pdf",sep=""),p3,height=4,width=7)
  
  p3b <- ggplot(AK_CA_Rec_Combo_meansOverSimReplicates[AK_CA_Rec_Combo_meansOverSimReplicates$label==dominance & AK_CA_Rec_Combo_meansOverSimReplicates$population=="AK",],aes(x=generation,y=meanLoad,color=population))+
    geom_line(data=AK_CA_PrePost[AK_CA_PrePost$label==dominance & AK_CA_PrePost$population=="AK",],aes(group=replicate,x=generation,y=L),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    geom_vline(data=hLines[hLines$population=="AK",],aes(xintercept = intercept),linetype="dashed",color="black")+
    facet_wrap(~population,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),legend.position = "none",strip.background = element_rect("transparent"))+
    xlab("generations since bottleneck")+
    ylab("mean genetic load")+
    scale_color_manual(values=unlist(colors))+
    ggtitle(paste(dominance," mutations",sep=""))
  
  p3b
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/RECOVERY.",dominance,".geneticLoad.AK.ONLY.pdf",sep=""),p3b,height=4,width=7)
  
}
##################### YOU ARE HERE --- NEED TO REDO LOAD CALC ########
################### III. 5 epoch double bottleneck *** NEED TO DO THIS CALCULATION FOR LOAD ####  #################
doubleBottleneckAK <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/poonehModels/5Epoch_20191204/AK/AK_5Epoch_data.txt",header=T)
doubleBottleneckAK$label <- NA
doubleBottleneckAK[doubleBottleneckAK$h==0,]$label <- "recessive"
doubleBottleneckAK[doubleBottleneckAK$h==0.5,]$label <- "additive"
# gather means:
doubleBottleneckAK_means = doubleBottleneckAK %>%
  group_by(model,population,label,generation)  %>%
  summarise(meanFitness=mean(W_meanFitness),meanLoad=mean(L))

################# IV. Model with Migration ###########
migModels <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/poonehModels/addingMigrationBetween_CA_AK/MigrationModels.data.txt",header=T)
AKbneckGen=8002 # after pops separated for 8K 
CAbneckGen=8012 # because CA bneck is 10 gens shorter
migModels_recent=migModels[(migModels$subpop==1 & migModels$generation>=AKbneckGen)|(migModels$subpop==2 & migModels$generation>=CAbneckGen),]
migModels_recent$genRescaled <- 0
migModels_recent[migModels_recent$subpop==2,]$genRescaled <- migModels_recent[migModels_recent$subpop==2,]$generation-CAbneckGen
migModels_recent[migModels_recent$subpop==1,]$genRescaled <- migModels_recent[migModels_recent$subpop==1,]$generation-AKbneckGen
##### get averages:
migModels_recent$label <- NA
migModels_recent[migModels_recent$h==0,]$label <- "recessive"
migModels_recent[migModels_recent$h==0.5,]$label <- "additive"
migModels_recent$pop <- NA
migModels_recent[migModels_recent$subpopulation==1,]$pop <- "AK"
migModels_recent[migModels_recent$subpopulation==2,]$pop <- "CA"
### nicer model names:
migModels_recent$modelID <- unlist(lapply(strsplit(as.character(migModels_recent$model),"\\."),tail,n=1)) # choose last element -- is element 4 for most models, but 3 for 'notranslocation'
# rename "notranslocation" to no migration
migModels_recent[migModels_recent$modelID=="NoTranslocation",]$modelID <- "NoGeneFlow"
#order factors:
migModels_recent$modelID <- factor(migModels_recent$modelID,levels=c("NoGeneFlow","1perGen","5perGen","10perGen","25for2Gen","25perGen"))
## the no translocation ones are bad 
# gather means:
migModels_recent_means = migModels_recent %>%
  group_by(model,pop,label,generation,modelID,genRescaled,subpopulation)  %>%
  summarise(meanFitness=mean(W_meanFitness),meanLoad=mean(L))
## when does migration start for each population?
hLines3=data.frame(pop=c("AK","AK","AK","CA","CA","CA"),intercept=c(0,35,54,0,25,43))

for(dominance in c("recessive","additive")){
  p5 <- ggplot(migModels_recent_means[migModels_recent_means$label==dominance & !(migModels_recent_means$modelID %in% c("5perGen","10perGen")),],aes(x=genRescaled,y=meanLoad,color=as.factor(modelID)))+
    #geom_line(data=migModels_recent[migModels_recent$label==dominance & !(migModels_recent$modelID %in% c("5perGen","10perGen")),],aes(group=interaction(rep,modelID),x=genRescaled,y=L_mutationLoad,color=as.factor(modelID)),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    geom_vline(data=hLines3,aes(xintercept = intercept),linetype="dashed",color="black")+
    facet_wrap(~pop,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),strip.background = element_rect("transparent"),legend.title = element_blank(),legend.position = "right",legend.background = element_rect(fill="transparent"),legend.text = element_text(size=rel(1.2)))+
    xlab("generations since population bottleneck")+
    ylab("mean genetic load")+
    ggtitle(paste(dominance," mutations",sep=""))+
    scale_color_manual(values=unlist(migrationColors))
  
  p5
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/MIGRATION.",dominance,".no5-10.geneticLoad.AK.CA.NoFineLines.pdf",sep=""),p5,height=3,width=8)
  
  p5b <- ggplot(migModels_recent_means[migModels_recent_means$label==dominance & !(migModels_recent_means$modelID %in% c("5perGen","10perGen")),],aes(x=genRescaled,y=meanLoad,color=as.factor(modelID)))+
    geom_line(data=migModels_recent[migModels_recent$label==dominance & !(migModels_recent$modelID %in% c("5perGen","10perGen")),],aes(group=interaction(rep,modelID),x=genRescaled,y=L,color=as.factor(modelID)),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    geom_vline(data=hLines3,aes(xintercept = intercept),linetype="dashed",color="black")+
    facet_wrap(~pop,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),strip.background = element_rect("transparent"),legend.title = element_blank(),legend.position = "right",legend.background = element_rect(fill="transparent"),legend.text = element_text(size=rel(1.2)))+
    xlab("generations since population bottleneck")+
    ylab("mean genetic load")+
    ggtitle(paste(dominance," mutations",sep=""))+
    scale_color_manual(values=unlist(migrationColors))
  
  p5b
  
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/MIGRATION.",dominance,".no5-10.geneticLoad.AK.CA.WithFineLines.pdf",sep=""),p5b,height=3,width=8)
  
  ###### plot just the alaska side: ##########
  p5c <- ggplot(migModels_recent_means[migModels_recent_means$label==dominance & !(migModels_recent_means$modelID %in% c("5perGen","10perGen")) & migModels_recent_means$pop=="AK",],aes(x=genRescaled,y=meanLoad,color=as.factor(modelID)))+
    #geom_line(data=migModels_recent[migModels_recent$label==dominance & !(migModels_recent$modelID %in% c("5perGen","10perGen")),],aes(group=interaction(rep,modelID),x=genRescaled,y=L_mutationLoad,color=as.factor(modelID)),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    geom_vline(data=hLines3[hLines3$pop=="AK",],aes(xintercept = intercept),linetype="dashed",color="black")+
    facet_wrap(~pop,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),strip.background = element_rect("transparent"),legend.title = element_blank(),legend.position = "right",legend.background = element_rect(fill="transparent"),legend.text = element_text(size=rel(1.2)))+
    xlab("generations since population bottleneck")+
    ylab("mean genetic load")+
    ggtitle(paste(dominance," mutations",sep=""))+
    scale_color_manual(values=unlist(migrationColors))
  
  p5c
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/MIGRATION.",dominance,".no5-10.geneticLoad.AK.Only.NoFineLines.pdf",sep=""),p5c,height=2.5,width=5)
  ########## can make same plot with fitness #########
  # p6 <- ggplot(migModels_recent_means[migModels_recent_means$label==dominance & !(migModels_recent_means$modelID %in% c("5perGen","10perGen")),],aes(x=genRescaled,y=meanFitness,color=as.factor(modelID)))+
  #   #geom_line(data=migModels_recent[migModels_recent$label==dominance & !(migModels_recent$modelID %in% c("5perGen","10perGen")),],aes(group=interaction(rep,modelID),x=genRescaled,y=L_mutationLoad,color=as.factor(modelID)),alpha=0.3,size=0.2)+
  #   geom_line(size=1.5)+
  #   theme_bw()+
  #   geom_vline(data=hLines3,aes(xintercept = intercept),linetype="dashed",color="black")+
  #   facet_wrap(~pop,ncol=2,scales="free_x")+
  #   theme(text = element_text(size=textsize),strip.background = element_rect("transparent"),legend.title = element_blank(),legend.position = "right",legend.background = element_rect(fill="transparent"),legend.text = element_text(size=rel(1.2)))+
  #   xlab("generations since population bottleneck")+
  #   ylab("mean fitness")+
  #   ggtitle(paste(dominance," mutations",sep=""))+
  #   scale_color_manual(values=unlist(migrationColors))
  # p6
  # ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/MIGRATION.",dominance,".no5-10.fitness.AK.CA.NoFineLines.pdf",sep=""),p6,height=3,width=8)
  # p6b <- ggplot(migModels_recent_means[migModels_recent_means$label==dominance & !(migModels_recent_means$modelID %in% c("5perGen","10perGen")),],aes(x=genRescaled,y=meanFitness,color=as.factor(modelID)))+
  #   geom_line(data=migModels_recent[migModels_recent$label==dominance & !(migModels_recent$modelID %in% c("5perGen","10perGen")),],aes(group=interaction(rep,modelID),x=genRescaled,y=W_meanFitness,color=as.factor(modelID)),alpha=0.3,size=0.2)+
  #   geom_line(size=1.5)+
  #   theme_bw()+
  #   geom_vline(data=hLines3,aes(xintercept = intercept),linetype="dashed",color="black")+
  #   facet_wrap(~pop,ncol=2,scales="free_x")+
  #   theme(text = element_text(size=textsize),strip.background = element_rect("transparent"),legend.title = element_blank(),legend.position = "right",legend.background = element_rect(fill="transparent"),legend.text = element_text(size=rel(1.2)))+
  #   xlab("generations since population bottleneck")+
  #   ylab("mean fitness")+
  #   ggtitle(paste(dominance," mutations",sep=""))+
  #   scale_color_manual(values=unlist(migrationColors))
  # p6b
  # ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/MIGRATION.",dominance,".no5-10.fitness.AK.CA.WithFineLines.pdf",sep=""),p6b,height=3,width=8)
  
  
}
############ Plot AK faceted by dominance:
migModels_recent_means$label2 <- paste(migModels_recent_means$label," mutations",sep="")
migModels_recent_means$label2 <- factor(migModels_recent_means$label2, levels=c("recessive mutations","additive mutations"))
p5d <- ggplot(migModels_recent_means[!(migModels_recent_means$modelID %in% c("5perGen","10perGen")) & migModels_recent_means$pop=="AK",],aes(x=genRescaled,y=meanLoad,color=as.factor(modelID)))+
  geom_line(size=1)+
  theme_bw()+
  geom_vline(data=hLines3[hLines3$pop=="AK",],aes(xintercept = intercept),linetype="dashed",color="black")+
  facet_wrap(~label2,ncol=1,scales="free_y")+
  theme(text = element_text(size=11),strip.background = element_blank(),legend.title = element_blank(),legend.position = "bottom",legend.justification = "left",legend.background = element_rect(fill="transparent"),legend.text = element_text(size=rel(0.9)),strip.text = element_text(hjust=0,size=textsize))+
  xlab("generations since bottleneck")+
  ylab("mean genetic load")+
  scale_color_manual(values=unlist(migrationColors))+
  guides(colour=guide_legend(ncol=2,nrow=2,byrow=TRUE))

p5d
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/MIGRATION.FacetDominance.no5-10.geneticLoad.AK.Only.NoFineLines.pdf",sep=""),p5d,height=4.8,width=3)

