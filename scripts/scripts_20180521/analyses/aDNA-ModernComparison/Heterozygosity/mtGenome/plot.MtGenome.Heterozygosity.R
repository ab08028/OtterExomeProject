########### Plotting heterozygosity on mt genome
# expectation: no hets
require(reshape2)
# two kinds of calls: hard calls where GP > = 0.95 or GP sums 

date="20191202-highcov-AFprior-MajorMinor4-mtGenomeOnly"
data.dir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/mtGenome/mtGenomeHeterozygosity/",date,"/",sep="")
GPInputs <- list.files(data.dir, pattern="mtGenomeOnly.txt",full.names = T)
HardCallInputs <- list.files(data.dir,pattern="HARDCALLSONLY.txt",full.names = T)

allGPInputs <- data.frame()
allHardCallInputs <- data.frame()
for(file in GPInputs){
  input <- read.table(file,header=T)
  allGPInputs <- rbind(allGPInputs,input)
}

for(file in HardCallInputs){
  input <- read.table(file,header=T)
  allHardCallInputs <- rbind(allHardCallInputs,input)
}

allGPInputs$minDepthLabel <- paste("min ",allGPInputs$Filter_PerIndividualDepthMinimum, " reads",sep="")
allGPInputs$minIndLabel <- paste("min ",allGPInputs$Filter_minIndsPerSite, " inds",sep="")

allHardCallInputs$minDepthLabel <- paste("min ",allHardCallInputs$Filter_PerIndividualDepthMinimum, " reads",sep="")
allHardCallInputs$minIndLabel <- paste("min ",allHardCallInputs$Filter_minIndsPerSite, " inds",sep="")

############ Plotting #########
require(ggplot2)

p1 <- ggplot(allGPInputs,aes(x=sample,y=callableSites,color=as.factor(Filter_ProbThresholdForCallableSite)))+
  geom_point()+
  facet_wrap(~minIndLabel~minDepthLabel)+
  coord_flip()+
  theme_bw()
p1
ggsave(paste(data.dir,"filters.callableSites.pdf",sep=""),p1,width=10,height=8)

p2a <- ggplot(allGPInputs,aes(x=sample,y=HetPerSite_TransversionsOnly,color=as.factor(Filter_ProbThresholdForCallableSite)))+
  geom_point()+
  facet_wrap(minIndLabel~minDepthLabel,ncol=3)+
  coord_flip()+
  theme_bw()
p2a
ggsave(paste(data.dir,"filters.hetPerSite.TVOnly.pdf",sep=""),p2a,width=10,height=8)

p2b <- ggplot(allGPInputs,aes(x=sample,y=HetPerSite,color=as.factor(Filter_ProbThresholdForCallableSite)))+
  geom_point()+
  facet_wrap(minIndLabel~minDepthLabel,ncol=3)+
  coord_flip()+
  theme_bw()
p2b

p3 <- ggplot(allGPInputs,aes(x=callableSites,y=HetPerSite_TransversionsOnly,color=as.factor(Filter_ProbThresholdForCallableSite)))+
  geom_point()+
  facet_wrap(minIndLabel~minDepthLabel,ncol=3,scales="free")+
  theme_bw()
p3
ggsave(paste(data.dir,"filters.callable.vs.het.TVOnly.pdf",sep=""),p3,width=10,height=8)

############# hard calls plots #########

p22 <- ggplot(allHardCallInputs,aes(x=sample,y=HetPerSite_TransversionsOnly))+
  geom_point()+
  facet_wrap(minIndLabel~minDepthLabel,ncol=3,scales="free_x")+
  coord_flip()+
  theme_bw()
p22

p32 <- ggplot(allHardCallInputs,aes(x=callableSites,y=HetPerSite_TransversionsOnly))+
  geom_point()+
  facet_wrap(minIndLabel~minDepthLabel,ncol=3,scales="free")+
  theme_bw()
p32
ggsave(paste(data.dir,"filters.callable.vs.het.TVOnly.HARDCALLs.pdf",sep=""),p32,width=10,height=8)


############## trying alternate ways to plot:
p22b <- ggplot(allHardCallInputs,aes(x=sample,y=HetPerSite_TransversionsOnly,fill=minDepthLabel))+
  geom_col(position="dodge")+
  facet_wrap(~minIndLabel,nrow=1,scales="free_x")+
  theme_bw()
p22b


p4 <- ggplot(allHardCallInputs,aes(x=sample,y=HetPerSite_TransversionsOnly,fill=minDepthLabel))+
  geom_col(position="dodge")+
  facet_wrap(~minIndLabel,nrow=1,scales="free_x")+
  theme_bw()
p4

# want to just do calls:
p5a <- ggplot(allHardCallInputs,aes(x=sample,y=sumHetGLsOrGPs_TransversionsOnly,fill=minDepthLabel))+
  geom_col(position="dodge")+
  facet_wrap(~minIndLabel,nrow=1,scales="free_x")+
  theme_bw()+
  ggtitle("Hard Calls")
p5a


p5b <- ggplot(allGPInputs,aes(x=sample,y=sumHetGLsOrGPs_TransversionsOnly,fill=minDepthLabel))+
  geom_col(position="dodge")+
  facet_wrap(~minIndLabel~Filter_ProbThresholdForCallableSite,nrow=1,scales="free_x")+
  theme_bw()+
  ggtitle("GP SUM")
p5b


p6 <- ggplot(allGPInputs,aes(x=sample,y=sumHetGLsOrGPs_TransversionsOnly,fill=minDepthLabel))+
  geom_col(position="dodge")+
  facet_wrap(~minIndLabel~Filter_ProbThresholdForCallableSite,nrow=1,scales="free_x")+
  theme_bw()+
  ggtitle("GP SUM")

###### let's plot them #######
bigInput <- read.table('/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/mtGenome/angsd-GLs/20191202-highcov-AFprior-MajorMinor4-mtGenomeOnly/angsdOut.mappedToelut.beagle.gprobs.gz',header=T)
countInput <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/aDNA-ModernComparison/mtGenome/angsd-GLs/20191202-highcov-AFprior-MajorMinor4-mtGenomeOnly/angsdOut.mappedToelut.counts.gz",header=T)
combo <- cbind(bigInput,countInput)
head(combo)
highHet <- combo[combo$Ind0.1 >=0.95 | combo$Ind1.1>=0.95 | combo$Ind2.1>=0.95,]
head(highHet)
ggplot(combo,aes(x=marker,y=Ind0.1))+
  geom_point()+
  geom_point(data=highHet,aes(x=marker,y=Ind1.1),color="red")+
  geom_point(data=highHet,aes(x=marker,y=Ind2.1),color="green")+
  geom_point(data=highHet,aes(x=marker,y=Ind3.1),color="blue")
require(reshape2)

# mark if transversion/transition;
# color by transitions/transversions
transversions=c("0,1",'1,0','0,3','3,0','1,2','2,1','2,3','3,2')
combo$alleles <- paste(combo$allele1,combo$allele2,sep=",")
combo$alleles %in% transversions
combo$TiTvLabel <- "NA"
combo[combo$alleles %in% transversions,]$TiTvLabel <- "Tv"
combo[!(combo$alleles %in% transversions),]$TiTvLabel <- "Ti"

comboMelt <- melt(combo)
justHets <- comboMelt[comboMelt$variable %in% c("Ind0.1","Ind1.1","Ind2.1","Ind3.1"),]
justHets$label <- "non-het"
justHets[justHets$value>=0.95,]$label <- "het"
ind0HetMarkers <- comboMelt[comboMelt$variable=="Ind0.1" & comboMelt$value>=0.95,]$marker
ind1HetMarkers <- comboMelt[comboMelt$variable=="Ind1.1" & comboMelt$value>=0.95,]$marker
ind2HetMarkers <- comboMelt[comboMelt$variable=="Ind2.1" & comboMelt$value>=0.95,]$marker
ind3HetMarkers <- comboMelt[comboMelt$variable=="Ind3.1" & comboMelt$value>=0.95,]$marker

p10 <- ggplot(justHets,aes(x=marker,y=value,color=label))+
  geom_col()+
  facet_wrap(~variable,ncol=1)+
  theme(axis.text.x = element_blank())
#p10
ggsave(paste(data.dir,"heterzyogistyAcrossMtGenome.pdf",sep=""),p10,width=10,height=8)
## want transversions!
justDepths <- comboMelt[comboMelt$variable %in% c("ind0TotDepth","ind1TotDepth","ind2TotDepth","ind3TotDepth"),]
justDepths$label <- "not het"

head(justDepths)
justDepths[justDepths$variable=="ind0TotDepth" & justDepths$marker %in% ind0HetMarkers,]$label <- "het"
justDepths[justDepths$variable=="ind1TotDepth" & justDepths$marker %in% ind1HetMarkers,]$label <- "het"
justDepths[justDepths$variable=="ind2TotDepth" & justDepths$marker %in% ind2HetMarkers,]$label <- "het"
justDepths[justDepths$variable=="ind3TotDepth" & justDepths$marker %in% ind3HetMarkers,]$label <- "het"
justDepths$individual <- as.factor(paste("Ind",unlist(lapply((strsplit(unlist(lapply(strsplit(as.character(justDepths$variable),split ="ind"),"[",2)),"Tot")),"[",1)),sep=""))
justHets$individual <- unlist(lapply(strsplit(as.character(justHets$variable),"\\."),"[",1))

p12 <- ggplot(justDepths,aes(x=marker,y=value))+
  geom_col()+
  geom_point(data=justHets[justHets$value>=0.95,],aes(x=marker,y=0,color=TiTvLabel),size=5)+
  facet_wrap(~individual,ncol=1)+
  theme(axis.text.x = element_blank())
#p10
ggsave(paste(data.dir,"depthAcrossMtGenome.TiTVsOnly.pdf",sep=""),p12,width=10,height=8)


# filter 
### get list of Tv Hets: ####
ind0TvHets <-combo[combo$Ind0.1 >= 0.95 & combo$TiTvLabel=="Tv",]

ind0TvHets$marker

# ref|NC_009692.1|_469   ref|NC_009692.1|_1304  ref|NC_009692.1|_7781  ref|NC_009692.1|_13387


############### Plot all hets ##########
allGPs <- comboMelt[!(comboMelt$variable %in% c("allele1","allele2","ind0TotDepth","ind1TotDepth","ind2TotDepth","ind3TotDepth")),]
allGPs$individual <- unlist(lapply(strsplit(as.character(allGPs$variable),"\\."),"[",1))
allGPs$variable2 <- "hom-ref"
allGPs[grep("\\.1",allGPs$variable),]$variable2 <- "het"
allGPs[grep("\\.2",allGPs$variable),]$variable2 <- "hom-alt"

p13 <- ggplot(allGPs[allGPs$value>=0.95 & allGPs$variable2!="hom-ref",],aes(x=marker,y=value,color=variable2,shape=TiTvLabel))+
  geom_point()+
  theme_bw()+
  facet_wrap(~individual,ncol=1)+
  theme(axis.text.x = element_blank())
p13
ggsave(paste(data.dir,"GPsAcrossGenome.pdf",sep=""),p13,width=10,height=8)
