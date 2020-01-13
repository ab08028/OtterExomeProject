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
AK_PrePost <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/poonehModels/3Epoch_20191204/AK/AK_data.txt",header=T)
AK_PrePost = AK_PrePost[AK_PrePost$generation<=36,]
### exclude recovery (doing that next )
CA_PrePost <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/poonehModels/3Epoch_20191204/CA/CA_data.txt",header=T)
CA_PrePost = CA_PrePost[CA_PrePost$generation<=26,]

AK_CA_PrePost_Combo <- rbind(AK_PrePost,CA_PrePost)
AK_CA_PrePost_Combo$label <- NA
AK_CA_PrePost_Combo[AK_CA_PrePost_Combo$h==0,]$label <- "recessive"
AK_CA_PrePost_Combo[AK_CA_PrePost_Combo$h==0.5,]$label <- "additive"
# gather means:
AK_CA_PrePost_Combo_meansOverSimReplicates = AK_CA_PrePost_Combo %>%
  group_by(model,population,label,generation)  %>%
  summarise(meanFitness=mean(W_meanFitness),meanLoad=mean(L_mutationLoad))

######################### p1: Plot the bottleneck load only ( no recovery ) ################
for(dominance in c("recessive","additive")){
  p1 <- ggplot(AK_CA_PrePost_Combo_meansOverSimReplicates[AK_CA_PrePost_Combo_meansOverSimReplicates$label==dominance,],aes(x=generation,y=meanLoad,color=population))+
    geom_line(data=AK_CA_PrePost_Combo[AK_CA_PrePost_Combo$label==dominance,],aes(group=rep,x=generation,y=L_mutationLoad),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    #geom_vline(data=hLines,aes(xintercept = intercept),linetype="dashed",color="darkorange")+
    facet_wrap(~population,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),legend.position = "none",strip.background = element_rect("transparent"))+
    xlab("generations since bottleneck")+
    ylab("mean genetic load")+
    scale_color_manual(values=unlist(colors))+
    ggtitle(paste(dominance," mutations",sep=""))
  p1
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/BOTTLENECK.",dominance,".geneticLoad.AK.CA.pdf",sep=""),p1,height=2,width=3)
  
  
  ########## p2: bottleneck fitness ###############
  # can make same plot with fitness 
  p2 <- ggplot(AK_CA_PrePost_Combo_meansOverSimReplicates[AK_CA_PrePost_Combo_meansOverSimReplicates$label==dominance,],aes(x=generation,y=meanFitness,color=population))+
    geom_line(data=AK_CA_PrePost_Combo[AK_CA_PrePost_Combo$label==dominance,],aes(group=rep,x=generation,y=W_meanFitness),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    #geom_vline(data=hLines,aes(xintercept = intercept),linetype="dashed",color="darkorange")+
    facet_wrap(~population,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),legend.position = "none",strip.background = element_rect("transparent"))+
    xlab("generations since bottleneck")+
    ylab("mean fitness")+
    scale_color_manual(values=unlist(colors))+
    ggtitle(paste(dominance," mutations",sep=""))
  p2
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/BOTTLENECK.",dominance,".fitness.AK.CA.pdf",sep=""),p2,height=2,width=3)
}

########## II. Pre-Post Contraction + Recovery #################
AK_Rec <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/poonehModels/3Epoch_20191204/AK/AK_data.txt",header=T)
### exclude recovery (doing that next )
CA_Rec <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/poonehModels/3Epoch_20191204/CA/CA_data.txt",header=T)

AK_CA_Rec_Combo <- rbind(AK_Rec,CA_Rec)
AK_CA_Rec_Combo$label <- NA
AK_CA_Rec_Combo[AK_CA_Rec_Combo$h==0,]$label <- "recessive"
AK_CA_Rec_Combo[AK_CA_Rec_Combo$h==0.5,]$label <- "additive"
# gather means:
AK_CA_Rec_Combo_meansOverSimReplicates = AK_CA_Rec_Combo %>%
  group_by(model,population,label,generation)  %>%
  summarise(meanFitness=mean(W_meanFitness),meanLoad=mean(L_mutationLoad))

######################### p3: Plot with recovery load ################
hLines=data.frame(population=c("AK","AK","CA","CA"),intercept=c(0,36,0,26))
for(dominance in c("recessive","additive")){
  
  p3 <- ggplot(AK_CA_Rec_Combo_meansOverSimReplicates[AK_CA_Rec_Combo_meansOverSimReplicates$label==dominance,],aes(x=generation,y=meanLoad,color=population))+
    geom_line(data=AK_CA_Rec_Combo[AK_CA_Rec_Combo$label==dominance,],aes(group=rep,x=generation,y=L_mutationLoad),alpha=0.3,size=0.2)+
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
  
  
  ########## p4: plot with recovery fitness #########
  p4 <- ggplot(AK_CA_Rec_Combo_meansOverSimReplicates[AK_CA_Rec_Combo_meansOverSimReplicates$label==dominance,],aes(x=generation,y=meanFitness,color=population))+
    geom_line(data=AK_CA_Rec_Combo[AK_CA_Rec_Combo$label==dominance,],aes(group=rep,x=generation,y=W_meanFitness),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    geom_vline(data=hLines,aes(xintercept = intercept),linetype="dashed",color="black")+
    facet_wrap(~population,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),legend.position = "none",strip.background = element_rect("transparent"))+
    xlab("generations since bottleneck")+
    ylab("mean fitness")+
    scale_color_manual(values=unlist(colors))+
    ggtitle(paste(dominance," mutations",sep=""))
  p4
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/RECOVERY.",dominance,".fitness.AK.CA.pdf",sep=""),p4,height=4,width=8)
}
################### III. 5 epoch double bottleneck #################
doubleBottleneckAK <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/poonehModels/5Epoch_20191204/AK/AK_5Epoch_data.txt",header=T)
doubleBottleneckAK$label <- NA
doubleBottleneckAK[doubleBottleneckAK$h==0,]$label <- "recessive"
doubleBottleneckAK[doubleBottleneckAK$h==0.5,]$label <- "additive"
# gather means:
doubleBottleneckAK_means = doubleBottleneckAK %>%
  group_by(model,population,label,generation)  %>%
  summarise(meanFitness=mean(W_meanFitness),meanLoad=mean(L_mutationLoad))

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
  summarise(meanFitness=mean(W_meanFitness),meanLoad=mean(L_mutationLoad))
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
    geom_line(data=migModels_recent[migModels_recent$label==dominance & !(migModels_recent$modelID %in% c("5perGen","10perGen")),],aes(group=interaction(rep,modelID),x=genRescaled,y=L_mutationLoad,color=as.factor(modelID)),alpha=0.3,size=0.2)+
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
  ########## can make same plot with fitness #########
  p6 <- ggplot(migModels_recent_means[migModels_recent_means$label==dominance & !(migModels_recent_means$modelID %in% c("5perGen","10perGen")),],aes(x=genRescaled,y=meanFitness,color=as.factor(modelID)))+
    #geom_line(data=migModels_recent[migModels_recent$label==dominance & !(migModels_recent$modelID %in% c("5perGen","10perGen")),],aes(group=interaction(rep,modelID),x=genRescaled,y=L_mutationLoad,color=as.factor(modelID)),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    geom_vline(data=hLines3,aes(xintercept = intercept),linetype="dashed",color="black")+
    facet_wrap(~pop,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),strip.background = element_rect("transparent"),legend.title = element_blank(),legend.position = "right",legend.background = element_rect(fill="transparent"),legend.text = element_text(size=rel(1.2)))+
    xlab("generations since population bottleneck")+
    ylab("mean fitness")+
    ggtitle(paste(dominance," mutations",sep=""))+
    scale_color_manual(values=unlist(migrationColors))
  p6
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/MIGRATION.",dominance,".no5-10.fitness.AK.CA.NoFineLines.pdf",sep=""),p6,height=3,width=8)
  p6b <- ggplot(migModels_recent_means[migModels_recent_means$label==dominance & !(migModels_recent_means$modelID %in% c("5perGen","10perGen")),],aes(x=genRescaled,y=meanFitness,color=as.factor(modelID)))+
    geom_line(data=migModels_recent[migModels_recent$label==dominance & !(migModels_recent$modelID %in% c("5perGen","10perGen")),],aes(group=interaction(rep,modelID),x=genRescaled,y=W_meanFitness,color=as.factor(modelID)),alpha=0.3,size=0.2)+
    geom_line(size=1.5)+
    theme_bw()+
    geom_vline(data=hLines3,aes(xintercept = intercept),linetype="dashed",color="black")+
    facet_wrap(~pop,ncol=2,scales="free_x")+
    theme(text = element_text(size=textsize),strip.background = element_rect("transparent"),legend.title = element_blank(),legend.position = "right",legend.background = element_rect(fill="transparent"),legend.text = element_text(size=rel(1.2)))+
    xlab("generations since population bottleneck")+
    ylab("mean fitness")+
    ggtitle(paste(dominance," mutations",sep=""))+
    scale_color_manual(values=unlist(migrationColors))
  p6b
  ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/manuscriptPlots/PrePostContraction_AK_CA/MIGRATION.",dominance,".no5-10.fitness.AK.CA.WithFineLines.pdf",sep=""),p6b,height=3,width=8)
  
  
}



