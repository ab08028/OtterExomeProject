require(ggplot2)
require(dplyr)
require(reshape2)
# Lets start with one population then can figure out how to generalize
############### fold SFS function ##########
# ####################### fold sfs function ###############
foldSFS <- function(sfs){
  foldedSFS <- data.frame()
  ss=length(sfs$frequency) - 1 # this is the ss in chromosomes
  foldedBin=ss/2  # half as many as ss ; do seq from 0 --> foldedBins to include monomorphics
  # note you have to start seq at 0 to include monomorphic bin
  for(i in seq(0,foldedBin)){
    # if it's not the center point (which doesn't get added together)
    # see wakeley coalescent eq 1.2
    if(i==ss-i){
      foldedSFS <- rbind.data.frame(foldedSFS, cbind.data.frame(frequency=i,count=(sfs[sfs$frequency==i,]$count + sfs[sfs$frequency==(ss-i),]$count)/2)) # if it's add midpoint (foldedLen, including 0 monomorphic bin), add together and divide by two (equivalent of not adding together)
    }
    else{
      foldedSFS <- rbind.data.frame(foldedSFS, cbind.data.frame(frequency=i,count=(sfs[sfs$frequency==i,]$count + sfs[sfs$frequency==(ss-i),]$count))) # if not at mid point, just add together like normal
    }
  }
  return(foldedSFS)
}

######## get empirical #######

pops=c("CA","AK","AL","COM","KUR")

todaysdate=format(Sys.Date(),format="%Y%m%d")
genotypeDate=20181119
projection.date=20181221
hetFilter=0.75 # level that hets were filtered
plot.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/SFS/",genotypeDate,"/easySFS_projection/projection-",projection.date,"/hetFilter-",hetFilter,"/",sep="")
dir.create(plot.dir,recursive = T)
syn.data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/cds/synonymous/projection-",projection.date,"-hetFilter-",hetFilter,"/fastsimcoal2/",sep="")
mis.data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/cds/missense/projection-",projection.date,"-hetFilter-",hetFilter,"/fastsimcoal2/",sep="")
neut.data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/",genotypeDate,"/easySFS_projection/neutral/projection-",projection.date,"-hetFilter-",hetFilter,"/fastsimcoal2-plusMonomorphic/",sep="")

all.empirical.sfs.folded.noMono <- data.frame()
for(pop in pops){
  
  ############## prepare the data ##############
  neut.input <- list.files(neut.data.dir,pattern=paste(pop,"_MAFpop0.obs",sep=""),full.names = T )
  syn.input <- list.files(syn.data.dir,pattern=paste(pop,"_MAFpop0.obs",sep=""),full.names = T )
  mis.input <- list.files(mis.data.dir,pattern=paste(pop,"_MAFpop0.obs",sep=""),full.names = T )
  neut.sfs <- melt(read.table(neut.input,skip = 1,header = T)) # skip first line: "13 folded "KUR""
  neut.sfs$mutLabel <- "neutral"
  syn.sfs <- melt(read.table(syn.input,skip = 1,header = T))
  syn.sfs$mutLabel <- "synonymous"
  mis.sfs <- melt(read.table(mis.input,skip = 1,header = T))
  mis.sfs$mutLabel <- "missense"
  # want to plot them separately and all together for the population
  # want to exclude d0_0 bin and 0-bins
  ########### put the syn, neut and mis sfses all together: ######
  
  all.sfs <- rbind.data.frame(syn.sfs,mis.sfs,neut.sfs)
  all.sfs$frequency <- sapply(strsplit(as.character(all.sfs$variable),"_"),"[",2)
  ######## exclude bins that are above the fold #########
  all.sfs_noMono_noZero <- all.sfs[all.sfs$value!=0 & all.sfs$frequency!="0",]
  ######## add proportional as well: ########
  all.sfs_noMono_noZero <- all.sfs_noMono_noZero %>%
    group_by(mutLabel) %>%
    mutate(proportion=value/sum(value))
  ####### order factors: #########
  all.sfs_noMono_noZero$mutLabel <- factor(all.sfs_noMono_noZero$mutLabel, levels=c("neutral","synonymous","missense"))
  ############ add identifiers to match simulated data #####
  all.sfs_noMono_noZero$population <- pop # population
  all.sfs_noMono_noZero$category <- "empirical"
  all.sfs_noMono_noZero$state <- "empirical"
  all.sfs_noMono_noZero$count <-   all.sfs_noMono_noZero$value
  
  ############ add to overall empirical set ########
  all.empirical.sfs.folded.noMono <- rbind(all.empirical.sfs.folded.noMono,data.frame(all.sfs_noMono_noZero))
  ########### plot count-sfses all together ########
  # sfsPlot <- ggplot(all.sfs_noMono_noZero,aes(x=as.numeric(frequency),y=proportion,fill=label))+
  #   geom_bar(stat="identity",position = "dodge")+
  #   theme_bw()+
  #   theme(legend.title=element_blank(),legend.background = element_rect(fill = "transparent"),legend.position = c(0.5,0.9),legend.direction = "horizontal",legend.spacing.x = unit(.15,"cm"),text=element_text(size=18))+
  #   xlab("frequency")+
  #   ggtitle(paste(pop," folded SFS (projected)\nhet filter: ",hetFilter,sep=""))+
  #   scale_x_continuous(breaks=c(seq(1,max(as.numeric(all.sfs_noMono_noZero$frequency)))))+
  #   scale_fill_manual(values=c("gray","dodgerblue","darkred"))
  # sfsPlot
  # ggsave(paste(plot.dir,pop,".cds.neutral.projected.prop.SFS.hetFilter.",hetFilter,".pdf",sep=""),sfsPlot,device = "pdf",height=5,width=7)
  # 
  # ########## Plot each individually ############
  # sfsPlot2 <- ggplot(all.sfs_noMono_noZero,aes(x=as.numeric(frequency),y=value,fill=label))+
  #   geom_bar(stat="identity",position = "dodge")+
  #   theme_bw()+
  #   theme(legend.title=element_blank(),legend.background = element_rect(fill = "transparent"),legend.position = "none",legend.direction = "horizontal",legend.spacing.x = unit(.15,"cm"),text=element_text(size=18))+
  #   xlab("frequency")+
  #   ylab("count")+
  #   ggtitle(paste(pop," folded SFS (projected)\nhet filter: ",hetFilter,sep=""))+
  #   scale_x_continuous(breaks=c(seq(1,max(as.numeric(all.sfs_noMono_noZero$frequency)))))+
  #   scale_fill_manual(values=c("gray","dodgerblue","darkred"))+
  #   facet_wrap(~label)
  # sfsPlot2
  # ggsave(paste(plot.dir,pop,".cds.neutral.projected.counts.faceted.SFS.hetFilter.",hetFilter,".pdf",sep=""),sfsPlot2,device = "pdf",height=5,width=12)
  # 
}
####### empirical sfses are now in: all.empirical.sfs.folded.noMono are folded and have no monomorphic counts ###########
################################## get simulated ###############################
sim.dir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/slim/cdsSimulations/SFSes/"
popmodels=c("COM/1D.3Epoch.1.5Mb.cds/20190404/","AK/1D.2Epoch.1.5Mb.cds/20190404/","AL/1D.2Epoch.1.5Mb.cds/20190404/") # skipping CA and KUR for now 
numreps=1 # number of reps 
# go through pre and post 
# empty df for count sfses
all.sim.sfs.folded.noMono <- data.frame()

for(pm in popmodels){ # go through each population-model dir
  for(state in c("Post","Pre")){  # do separtely for pre and post contraction 
    for(rep in seq(1,numreps)){ # go through each replicate
      for(mutType in c(1,2)){ # go through each mutation type
    # mutation type 1 
    sim.file <- list.files(paste(sim.dir,pm,sep=""),pattern=paste("*rep.",rep,".",state,"Contraction.slim.mutType",as.character(mutType),".unfolded.sfs.R.format",sep=""),full.names = T) # this is for one rep, but could potentially do all reps and average here by turning rep into *
    # Fill in info: 
    ################## mutation type 1 ####################
    sfs <- read.table(sim.file,header=T)
    sfs.folded <- foldSFS(sfs) #### the input SFS from slim is UNFOLDED, fold here ***
    # get rid of monomorphic sites
    sfs.folded.noMono <- sfs.folded[sfs.folded$frequency!=0,]
    sfs.folded.noMono
    sfs.folded.noMono$mutType <- mutType
    sfs.folded.noMono$popModel <- pm
    sfs.folded.noMono$population <- unlist(lapply(strsplit(sfs.folded.noMono$popModel,"/"),"[",1))
    sfs.folded.noMono$model <- unlist(lapply(strsplit(sfs.folded.noMono$popModel,"/"),"[",2))
    sfs.folded.noMono$modelDate <- unlist(lapply(strsplit(sfs.folded.noMono$popModel,"/"),"[",3))
    sfs.folded.noMono$state <- paste("simulated: ",state,"-Contraction",sep="")
    sfs.folded.noMono$category <- "simulated"
    # also get proportions: 
    sfs.folded.noMono <- sfs.folded.noMono %>%
      mutate(proportion=count/sum(count))
    all.sim.sfs.folded.noMono <- rbind(all.sim.sfs.folded.noMono, sfs.folded.noMono)
    }
  }
}
}

####### simulated sfses are now in: all.sim.sfs.folded.noMono and are folded and have no monomorphic counts ###########
## Based on my slim script, mutation type 1 is "synonymous" and mutation type 2 is "missense" (not this is different from Chris' script labels)
# so assign those labels
all.sim.sfs.folded.noMono$mutLabel <- NA
all.sim.sfs.folded.noMono[all.sim.sfs.folded.noMono$mutType=="1",]$mutLabel <- "synonymous"
all.sim.sfs.folded.noMono[all.sim.sfs.folded.noMono$mutType=="2",]$mutLabel <- "missense"
# skipping neutral for now, but should find those simulated neutral sfses

################## combine empirical and simulated cds #################
# things in common that I need to plot are:
# "mutLabel"   "frequency"  "proportion" "pop"        "category"   "state"      "count"   
colsIWant <- intersect(names(all.empirical.sfs.folded.noMono),names(all.sim.sfs.folded.noMono))
# and want to exclude some pops:
popsIWant <- c("AK","AL") # skip COM for now -- it's weird

combo.cds <- rbind(all.empirical.sfs.folded.noMono[all.empirical.sfs.folded.noMono$population %in% popsIWant,colsIWant],all.sim.sfs.folded.noMono[all.sim.sfs.folded.noMono$population %in% popsIWant,colsIWant])
# exclude neutral for now --> 
# order factors
####### order factors: #########
combo.cds$state <- factor(combo.cds$state, levels=c("empirical","simulated: Pre-Contraction","simulated: Post-Contraction"))

p1 <- ggplot(combo.cds[combo.cds$mutLabel!="neutral",],aes(x=as.numeric(frequency),y=proportion,fill=state))+
  geom_bar(position='dodge',stat="identity")+
  facet_grid(mutLabel~as.factor(population),scales="free")+
  theme_bw()+
  scale_fill_manual(values=c("darkgray","blue","darkred"))+
  xlab("Frequency")+
  ylab("Proportion")+
  scale_x_continuous(breaks=c(seq(1,max(as.numeric(combo.cds$frequency)))))
p1
# makea  better output dir: 
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/compareSimulatedEmpiricalSFSes/test.ComparingSimulatedEmpiricalCDS.SFSes.AK.AL.COM.",todaysdate,".pdf",sep=""),p1,device="pdf",width=11,height=5)


# plot ratio of NS and S 
NS_S_Ratio <- combo.cds %>%
  group_by(state,population,mutLabel) %>%
  tally(sum(count)) %>%
  summarise(Ratio = n[mutLabel == "missense"] / n[mutLabel == "synonymous"])
ggplot(NS_S_Ratio,aes(x=population,y=Ratio,fill=state))+
  geom_bar(stat="identity", position="dodge")+
  ylab("Ratio of Missense / Synonymous SNPs")+
  theme_bw()+
  scale_fill_manual(values=c("darkgray","blue","darkred"))
# bad DFE? 